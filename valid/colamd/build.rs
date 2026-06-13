fn main() {
    println!("cargo::rerun-if-changed=ref/colamd.h");
    println!("cargo::rerun-if-changed=ref/colamd.c");
    cc::Build::new().include("ref").file("ref/colamd.c").flags(["-std=c99", "-O2", "-Wall", "-Wextra"]).compile("reference");
}
