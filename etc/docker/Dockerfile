# syntax=docker/dockerfile:1

ARG BASE_IMAGE_VERSION=trixie

FROM debian:${BASE_IMAGE_VERSION} AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG MAKE_JOBS=4
ARG TARGETARCH

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        autoconf \
        automake \
        autotools-dev \
        bison \
        ca-certificates \
        cmake \
        file \
        flex \
        g++ \
        gcc \
        gfortran \
        git \
        libfftw3-dev \
        libglut-dev \
        libhdf5-dev \
        libopenblas-dev \
        m4 \
        make \
        patch \
        perl \
        pkg-config \
        python3 \
        unzip \
        wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

# These legacy bundled packages do not build cleanly on ARM64. PETSc still
# provides the current MMG and ParMMG libraries on that architecture.
RUN set -eux; \
    set --; \
    if [ "${TARGETARCH:-}" = "arm64" ]; then \
        set -- --disable-mmg3d --disable-nlopt --disable-ipopt; \
    fi; \
    autoreconf -i; \
    ./configure \
        --enable-download \
        --enable-optim \
        --enable-generic \
        --prefix=/usr/local \
        "$@"

RUN ./3rdparty/getall -a -o PETSc \
    && make -C 3rdparty/ff-petsc petsc-slepc \
    && ./reconfigure

RUN make -j"${MAKE_JOBS}" \
    && make install \
    && find /usr/local/lib /usr/local/ff-petsc \
        -type f \( -name '*.a' -o -name '*.la' \) \
        -delete

FROM debian:${BASE_IMAGE_VERSION}-slim AS runtime

ARG DEBIAN_FRONTEND=noninteractive
ARG FREEFEM_VERSION=unknown

LABEL org.opencontainers.image.title="FreeFEM" \
      org.opencontainers.image.description="FreeFEM finite element software" \
      org.opencontainers.image.url="https://freefem.org" \
      org.opencontainers.image.source="https://github.com/FreeFem/FreeFem-sources" \
      org.opencontainers.image.version="${FREEFEM_VERSION}"

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        ca-certificates \
        g++ \
        gcc \
        gfortran \
        libfftw3-double3 \
        libglut3.12 \
        libhdf5-310 \
        libopenblas0-pthread \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/ /usr/local/

ENV PATH="/usr/local/bin:${PATH}" \
    GLOB="*.edp"

WORKDIR /data

RUN printf '%s\n' 'cout << "FreeFEM container smoke test" << endl;' > /tmp/smoke.edp \
    && FreeFem++-nw -v 0 /tmp/smoke.edp \
    && rm /tmp/smoke.edp

CMD ["bash", "-c", "shopt -s nullglob; files=(/data/$GLOB); for file in \"${files[@]}\"; do ff-mpirun -n 1 \"$file\" || exit $?; done"]
