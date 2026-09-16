# syntax=docker/dockerfile:1
FROM public.ecr.aws/docker/library/ubuntu:26.04 AS build

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_JOBS=2
ENV VCPKG_ROOT=/opt/vcpkg
ENV VCPKG_INSTALLED_DIR=/opt/vcpkg/installed
ENV _VCPKG_INSTALLED_DIR=/opt/vcpkg/installed
ENV VCPKG_MAX_CONCURRENCY=${BUILD_JOBS}

# These packages bootstrap vcpkg and provide the host tools required by gmsh[occ].
RUN apt-get update && apt-get install --yes --no-install-recommends \
    autoconf \
    autoconf-archive \
    automake \
    build-essential \
    ca-certificates \
    cmake \
    curl \
    git \
    libalglib-dev \
    libann-dev \
    libblas-dev \
    libcgns-dev \
    libfontconfig1-dev \
    libfltk1.3-dev \
    libfreetype-dev \
    libgmsh-dev \
    libgl2ps-dev \
    libglu1-mesa-dev \
    libglut-dev \
    libgmp-dev \
    libjpeg-dev \
    libtool \
    libx11-dev \
    libxext-dev \
    libxi-dev \
    libxmu-dev \
    libxt-dev \
    libmetis-dev \
    libocct-data-exchange-dev \
    libocct-foundation-dev \
    libocct-modeling-algorithms-dev \
    libocct-modeling-data-dev \
    libocct-ocaf-dev \
    libpng-dev \
    mesa-common-dev \
    ninja-build \
    patchelf \
    pkg-config \
    python3 \
    tar \
    unzip \
    zip \
    && ln -s "$(ldconfig -p | perl -ne 'print $1 and exit if /libvoro\+\+\.so.* => (.*)/')" /usr/lib/x86_64-linux-gnu/libvoro++.so \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/microsoft/vcpkg.git "$VCPKG_ROOT" \
    && git -C "$VCPKG_ROOT" checkout eed289f6e06a5e7a5c9e6a729671b0a56af7dd69 \
    && "$VCPKG_ROOT/bootstrap-vcpkg.sh" -disableMetrics

WORKDIR /src
COPY vcpkg.json ./
COPY triplets ./triplets
RUN "$VCPKG_ROOT/vcpkg" install --triplet x64-linux --overlay-triplets=/src/triplets

COPY . .

RUN cmake --preset gnu -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build --target tulip --parallel "$BUILD_JOBS"

# Package the executable with every shared library resolved by the dynamic linker.
RUN mkdir -p /bundle/bin /bundle/lib \
    && cp build/bin/tulip /bundle/bin/ \
    && ldd build/bin/tulip | awk '/=> \/[^ ]+/ { print $3 } /^\// { print $1 }' | sort -u | xargs -r -I{} cp {} /bundle/lib/ \
    && patchelf --set-rpath '$ORIGIN/../lib' /bundle/bin/tulip \
    && for library in /bundle/lib/*; do patchelf --set-rpath '$ORIGIN' "$library" || true; done

FROM scratch AS artifact
COPY --from=build /bundle/ /
