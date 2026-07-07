ARG IMAGE=intersystems/iris-community:latest-em
FROM $IMAGE

WORKDIR /home/irisowner/dev
COPY . .

## Python dep: Polars (vectorized parse / shared tail). No NumPy — the C kernel
## writes into plain ctypes buffers that Polars reads via memoryview.
RUN pip install --no-cache-dir --break-system-packages polars

## Build the C scan kernel (libdeflate gunzip + single-pass min/max, OpenMP).
## Best-effort: if gcc/libdeflate are missing the build simply skips the .so and
## Gaia.Benchmark falls back to the pure-Polars engine at run time. Runs as root
## so it can write into src/ (owned by root after COPY), then hands the .so back
## to irisowner. The volume mount shadows this at run time, so gaia_fast also
## lazily recompiles if the .so isn't present.
USER root
RUN LIBDEFLATE="$(ls /usr/lib/*/libdeflate.so.0 2>/dev/null | head -1)"; \
    if command -v gcc >/dev/null && [ -n "$LIBDEFLATE" ]; then \
        { gcc -O3 -march=native -funroll-loops -fopenmp -fPIC -shared src/gaiascan.c "$LIBDEFLATE" -o src/gaiascan.so || \
          gcc -O3 -fopenmp -fPIC -shared src/gaiascan.c "$LIBDEFLATE" -o src/gaiascan.so; } && \
        chown irisowner:irisowner src/gaiascan.so && \
        echo "built gaiascan.so"; \
    else \
        echo "gcc/libdeflate unavailable — skipping C kernel (Polars fallback)"; \
    fi
USER irisowner

## Embedded Python environment
ENV IRISUSERNAME="_SYSTEM"
ENV IRISPASSWORD="SYS"
ENV IRISNAMESPACE="USER"
ENV PYTHON_PATH=/usr/irissys/bin/
ENV PATH="/usr/irissys/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/irisowner/bin"

RUN --mount=type=bind,src=.,dst=. \
    iris start IRIS && \
	iris merge IRIS merge.cpf && \
	iris session IRIS < iris.script && \
    iris stop IRIS quietly safely