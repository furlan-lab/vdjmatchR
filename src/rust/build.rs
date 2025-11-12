fn main() {
    // Generate R wrappers at build time, locating the package root by DESCRIPTION.
    extendr_build::configure();
}
