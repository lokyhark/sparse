set -ex
cargo check --workspace --all-features --all-targets --release
cargo +nightly fmt --all -- --check
cargo clippy --workspace --all-features --all-targets --release -- --deny warnings
cargo build --workspace --all-features --all-targets --release
cargo test --workspace --all-features --all-targets --release
cargo +nightly miri test --workspace --all-features --all-targets
cargo doc --workspace --all-features --no-deps --release
cargo test --workspace --all-features --doc
