use std::{
    error::Error,
    ffi::OsStr,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

use flate2::read::GzDecoder;
use tar::Archive;

unsafe extern "C" {
    pub fn colamd_recommended(nnz: i32, nrows: i32, ncols: i32) -> i32;
    pub fn colamd(nrows: i32, ncols: i32, len: i32, a: *mut i32, p: *mut i32, knobs: *mut f64, stats: *mut i32) -> i32;
}

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::var("SUITE_SPARSE_COLLECTION_PATH").expect("SUITE_SPARSE_COLLECTION_PATH variable not defined");
    let path = PathBuf::from(path);
    let mut paths = Vec::new();
    visit(&path, &mut paths)?;
    println!("| {:^40} | {:^8} | {:4} | {:^20} | {:^20} | {:^20} | {:^8} |", "name", "id", "kind", "nrows", "ncols", "nnz", "valid");
    println!("| {0:-^40} | {0:-^8} | {0:-<4} | {0:-^20} | {0:-^20} | {0:-^20} | {0:-^8} |", "");
    for path in paths {
        match path.extension() {
            Some(ext) if ext == OsStr::new("gz") => (),
            _ => continue,
        }

        let file = File::open(path)?;
        let tar = GzDecoder::new(file);
        let mut archive = Archive::new(tar);
        for entry in archive.entries()? {
            let entry = match entry {
                Ok(entry) => entry,
                Err(_) => continue,
            };
            let path = match entry.path() {
                Ok(path) => path,
                Err(_) => continue,
            };
            match path.extension() {
                Some(ext) if ext == OsStr::new("rb") => (),
                _ => continue,
            }
            let mut buf = BufReader::new(entry);
            let mut line = String::new();
            buf.read_line(&mut line)?;
            let title = match line.get(..72) {
                Some(slice) => match slice.split(';').next() {
                    Some(slice) => slice.trim().trim_end_matches('|').trim().to_owned(),
                    None => panic!("title not found"),
                },
                None => panic!("title not found"),
            };
            let key = match line.get(72..) {
                Some(slice) => slice.trim().trim_start_matches('|').to_owned(),
                None => panic!("key not found"),
            };
            line.clear();
            buf.read_line(&mut line)?;
            let lptr: usize = match line.get(14..28) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse lptr")),
                None => panic!("lptr not found"),
            };
            let lind: usize = match line.get(28..42) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse lind")),
                None => panic!("lind not found"),
            };
            line.clear();
            buf.read_line(&mut line)?;
            let kind = match line.get(..3) {
                Some(slice) => slice,
                None => panic!("kind not found"),
            };

            let nrows: i32 = match line.get(14..28) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse nrows")),
                None => panic!("nrows not found"),
            };
            let ncols: i32 = match line.get(28..42) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse ncols")),
                None => panic!("ncols not found"),
            };
            let nnz: i32 = match line.get(42..56) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse nnz")),
                None => panic!("nnz not found"),
            };
            print!("| {title:40} | {key:>8} | {kind:4} | {nrows:20} | {ncols:20} | {nnz:20} |");
            std::io::stdout().flush().unwrap();
            match kind.get(1..3) {
                Some("ua") => (),
                _ => {
                    println!(" {:^7} |", "🟣");
                    continue;
                }
            }
            buf.read_line(&mut line)?;
            line.clear();
            for _ in 0..(lptr) {
                buf.read_line(&mut line)?;
            }
            let mut p = Vec::new();
            for ptr in line.split_whitespace() {
                p.push(ptr.parse::<i32>().unwrap() - 1);
            }

            let len = unsafe { colamd_recommended(nnz, nrows, ncols) };
            if len == 0 {
                panic!("colamd recommended error");
            }
            let mut a = Vec::with_capacity(len as usize);
            line.clear();
            for _ in 0..(lind) {
                buf.read_line(&mut line)?;
            }
            for ind in line.split_whitespace() {
                a.push(ind.parse::<i32>().unwrap() - 1);
            }

            // Implementation
            let order = sparse::colamd::colamd(nrows, ncols, &p, &a[..nnz as usize])?;
            // Reference
            let mut knobs = vec![0.0; 20];
            knobs[0] = 10.0; // Dense row control
            knobs[1] = 10.0; // Dense col control
            knobs[2] = 1.0; // Agressive
            let mut stats = vec![0; 20];
            let ok = unsafe { colamd(nrows, ncols, len, a.as_mut_ptr(), p.as_mut_ptr(), knobs.as_mut_ptr(), stats.as_mut_ptr()) };
            if ok == 0 {
                panic!("colamd error");
            }

            // // Check
            let mut valid = true;
            for i in 0..ncols as usize {
                if p[i] != order.order()[i] {
                    valid = false;
                    println!("invalid order at index {i}: ref = {}, impl = {}", p[i], order.order()[i]);
                    break;
                }
            }
            if valid {
                println!(" {:^7} |", "🟢");
            } else {
                println!(" {:^7} |", "🔴");
            }
        }
    }
    Ok(())
}

fn visit(directory: &Path, collection: &mut Vec<PathBuf>) -> std::io::Result<()> {
    if directory.is_file() {
        collection.push(directory.to_owned());
    } else if directory.is_dir() {
        for entry in std::fs::read_dir(directory)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                visit(&path, collection)?;
            } else {
                collection.push(path);
            }
        }
    } else {
        panic!("invalid path {}", directory.display());
    }
    Ok(())
}
