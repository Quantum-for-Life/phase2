extern crate cmake;

fn main() {
    let mut config = cmake::Config::new("../../");

    config.define("DISTRIBUTED", "OFF");

    #[cfg(feature = "mpi")]
    config.define("DISTRIBUTED", "ON");

    // This will also compile QuEST and log
    config.build_target("phase2");
    let dst = config.build();
    println!(
        "cargo:rustc-link-search=native={}/build/lib/QuEST/QuEST",
        dst.display()
    );
    println!(
        "cargo:rustc-link-search=native={}/build/src/phase2",
        dst.display()
    );
    println!("cargo:rustc-link-lib=dylib=phase2");
}
