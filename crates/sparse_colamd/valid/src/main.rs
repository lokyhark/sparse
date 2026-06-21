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
    // Retrieve command line arguments.
    let mut args = std::env::args();
    // Ignore path of the executable.
    args.next();
    // Expected Arguments.
    let mut path = None; // Path to sparse matrice archives (required).
    let mut max_bytes = 2u64.pow(30); // Maximum file size (default = 1GiB).
    let mut max_nnz = 10_000_000; // Maximum number of non zeros.
    let mut max_rows = 10_000_000; // Maximum number of rows.
    let mut max_cols = 10_000_000; // Maximum number of cols.
    // Parse command line arguments.
    while let Some(arg) = args.next() {
        match arg.as_str() {
            // Sparse matrice archives path.
            "-p" | "--path" => match args.next() {
                Some(arg) => {
                    let arg = PathBuf::from(arg);
                    if arg.exists() {
                        path = Some(arg)
                    } else {
                        return Err("invalid sparse matrice archives path".into());
                    }
                }
                None => return Err("sparse matrice archives path not given".into()),
            },
            // Maximum number of bytes to read.
            "-b" | "--bytes" => match args.next() {
                Some(arg) => {
                    max_bytes = match arg.parse() {
                        Ok(bytes) => bytes,
                        Err(error) => return Err(format!("failed to parse maximum bytes: {error}").into()),
                    }
                }
                None => return Err("maximum bytes size not given".into()),
            },
            // Maximum number of non zeros to consider.
            "-n" | "--nnz" => match args.next() {
                Some(arg) => {
                    max_nnz = match arg.parse() {
                        Ok(nnz) => nnz,
                        Err(error) => return Err(format!("failed to parse number of non zeros: {error}").into()),
                    }
                }
                None => return Err("maximum number of non zero not given".into()),
            },
            // Maximum number of rows to consider.
            "-r" | "--rows" => match args.next() {
                Some(arg) => {
                    max_rows = match arg.parse() {
                        Ok(rows) => rows,
                        Err(error) => return Err(format!("failed to parse number of rows: {error}").into()),
                    }
                }
                None => return Err("maximum number of rows not given".into()),
            },
            // Maximum number of cols to consider.
            "-c" | "--cols" => match args.next() {
                Some(arg) => {
                    max_cols = match arg.parse() {
                        Ok(cols) => cols,
                        Err(error) => return Err(format!("failed to parse number of cols: {error}").into()),
                    }
                }
                None => return Err("maximum number of cols not given".into()),
            },
            // Invalid command line arguments.
            _ => return Err(format!("invalid command line argument: {arg}").into()),
        }
    }

    // Check that path is set.
    let path = match path {
        Some(path) => path,
        None => return Err("sparse matrice archives path not defined".into()),
    };
    // Walk sparse matrice archives.
    let mut paths = Vec::new();
    visit(&path, &mut paths)?;

    // Print table header.
    println!("| {:^40} | {:^8} | {:4} | {:^20} | {:^20} | {:^20} | {:^8} |", "name", "id", "kind", "nrows", "ncols", "nnz", "valid");
    println!("| {0:-^40} | {0:-^8} | {0:-<4} | {0:-^20} | {0:-^20} | {0:-^20} | {0:-^8} |", "");
    // Iterate over all archives.
    for path in paths {
        // Ignore irrelevant files.
        match path.extension() {
            Some(ext) if ext == OsStr::new("gz") => (),
            _ => continue,
        }
        // Check file size.
        // This check was added to prevent loading extremely large matrices available in the suite sparse matrix collection.
        let file = File::open(&path)?;
        let metadata = file.metadata()?;
        if metadata.len() > max_bytes {
            let name = path.file_name().expect("failed to retrieve file name").to_str().expect("failed to parse file name");
            println!("| {1:<40} | {0:^8} | {0:<4} | {0:^20} | {0:^20} | {0:^20} | {2:^7} |", "", name, "🟣");
            continue;
        }
        // Load archive.
        let tar = GzDecoder::new(file);
        let mut archive = Archive::new(tar);
        // Iterate over archive entries.
        for entry in archive.entries()? {
            // Retrieve full entry path.
            let entry = match entry {
                Ok(entry) => entry,
                Err(_) => continue,
            };
            let path = entry.path()?;
            // Consider only rutherford boeing format.
            match path.extension() {
                Some(ext) if ext == OsStr::new("rb") => (),
                _ => continue,
            }
            // Bufread to increase performances.
            let mut buf = BufReader::new(entry);
            // Allocate line.
            let mut line = String::new();
            // Parse matrice generale informations.
            buf.read_line(&mut line)?;
            let title = match line.get(..72) {
                Some(slice) => match slice.split(';').next() {
                    Some(slice) => slice.trim().trim_end_matches('|').trim().to_owned(),
                    None => return Err("title not found".into()),
                },
                None => return Err("title not found".into()),
            };
            print!("| {title:40} ");
            let key = match line.get(72..) {
                Some(slice) => slice.trim().trim_start_matches('|').to_owned(),
                None => return Err("key not found".into()),
            };
            line.clear();
            buf.read_line(&mut line)?;
            let lptr: usize = match line.get(14..28) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse lptr")),
                None => return Err("lptr not found".into()),
            };
            let lind: usize = match line.get(28..42) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse lind")),
                None => return Err("lind not found".into()),
            };
            line.clear();
            buf.read_line(&mut line)?;
            let kind = match line.get(..3) {
                Some(slice) => slice,
                None => return Err("kind not found".into()),
            };

            let nrows: i32 = match line.get(14..28) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse nrows")),
                None => return Err("nrows not found".into()),
            };
            let ncols: i32 = match line.get(28..42) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse ncols")),
                None => return Err("ncols not found".into()),
            };
            let nnz: i32 = match line.get(42..56) {
                Some(slice) => slice.trim().parse().unwrap_or_else(|_| panic!("failed to parse nnz")),
                None => return Err("nnz not found".into()),
            };
            // Print matrix information
            print!("| {key:>8} | {kind:4} | {nrows:20} | {ncols:20} | {nnz:20} |");
            // Flush stdout to force display before processing.
            std::io::stdout().flush().unwrap();
            // Check matrix size.
            if nnz > max_nnz || nrows > max_rows || ncols > max_cols {
                println!(" {:^7} |", "🟣");
                continue;
            }
            // Ignore FORTRAN format lines.
            buf.read_line(&mut line)?;
            // Parse column pointers.
            line.clear();
            for _ in 0..(lptr) {
                buf.read_line(&mut line)?;
            }
            let mut p = Vec::new();
            for ptr in line.split_whitespace() {
                match ptr.parse::<i32>() {
                    // Pointers are 1-based in RB format.
                    Ok(ptr) => p.push(ptr - 1),
                    Err(error) => return Err(format!("failed to parse ptr: {error}").into()),
                }
            }

            // Determine indices workspace size according to reference colamd.
            let len = unsafe { colamd_recommended(nnz, nrows, ncols) };
            if len == 0 {
                return Err("colamd recommended error".into());
            }
            // Allocate indices workspace.
            let mut a = Vec::with_capacity(len as usize);
            line.clear();
            // Parse row indices.
            for _ in 0..(lind) {
                buf.read_line(&mut line)?;
            }
            for ind in line.split_whitespace() {
                match ind.parse::<i32>() {
                    // Indices are 1-based in RB format.
                    Ok(ind) => a.push(ind - 1),
                    Err(error) => return Err(format!("failed to parse index: {error}").into()),
                }
            }

            // Implementation calculation.
            let order = match sparse::colamd::colamd(nrows, ncols, &p, &a[..nnz as usize]) {
                Ok(order) => order,
                Err(_) => {
                    println!(" {:^7} |", "🔴");
                    continue;
                }
            };

            // Reference calculation.
            let mut knobs = vec![0.0; 20];
            knobs[0] = 10.0; // Dense row control
            knobs[1] = 10.0; // Dense col control
            knobs[2] = 1.0; // Agressive
            let mut stats = vec![0; 20];
            let ok = unsafe { colamd(nrows, ncols, len, a.as_mut_ptr(), p.as_mut_ptr(), knobs.as_mut_ptr(), stats.as_mut_ptr()) };
            if ok == 0 {
                return Err("colamd error".into());
            }

            // Check
            let mut valid = true;
            for (i, pos) in p.iter().enumerate().take(ncols as usize) {
                if *pos != order.order()[i] {
                    valid = false;
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
