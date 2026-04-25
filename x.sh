set -ex
cargo check --workspace --all-features --all-targets --release
cargo fmt --all -- --check
cargo clippy --workspace --all-features --all-targets --release -- --deny warnings
cargo build --workspace --all-features --all-targets --release
cargo test --workspace --all-features --all-targets --release
cargo doc --workspace --all-features --no-deps --release
cargo test --workspace --all-features --doc
