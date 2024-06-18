// We need to forward routine registration from C to Rust
// to avoid the linker removing the static library.

void R_init_apt_rust_extendr(void *dll);

void R_init_apt_rust(void *dll) {
    R_init_apt_rust_extendr(dll);
}
