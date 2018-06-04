extern crate cmake;

use cmake::Config;
use std::env;

fn main() {
    let dst = Config::new("libdivsufsort")
        .define("BUILD_EXAMPLES", "OFF")
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("BUILD_DIVSUFSORT64", "ON")
        .define("BUILD_SHARED_LIBS", "OFF")
        .define("CMAKE_INSTALL_LIBDIR", env::var("OUT_DIR").unwrap())
        .build();
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=static=divsufsort64");
}
