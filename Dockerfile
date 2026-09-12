
# image for rust
FROM rust:latest as builder

WORKDIR /app

COPY . .

RUN cargo build --release

FROM debian:bookworm-slim

# runtime
RUN apt-get update && apt-get install -y \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/vcf-reformatter /usr/local/bin/vcf-reformatter

ENTRYPOINT ["vcf-reformatter"]