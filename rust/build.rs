extern crate cmake;

fn main() {

    let mut config = cmake::Config::new("..");
    // This will also compile QuEST
    config.build_target("phase2");
    let dst = config.build();
    println!(
        "cargo:rustc-link-search=native={}/build/lib/QuEST/QuEST",
        dst.display()
    );
    println!(
        "cargo:rustc-link-search=native={}/build/src/",
        dst.display()
    );
    println!("cargo:rustc-link-lib=dylib=phase2");
}
