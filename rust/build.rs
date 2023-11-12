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
//     println!("cargo:rustc-link-lib=static=QuEST");
    println!("cargo:rustc-link-lib=dylib=circ");
    println!("cargo:rustc-link-lib=dylib=linen");
    println!("cargo:rustc-link-lib=dylib=rayon");
}
