#!/bin/bash
# start.sh — сборка и запуск AAMCC (GRATE)
#
# Использование:
#   ./start.sh            # собрать и запустить run_afterburner.C (по умолчанию)
#   ./start.sh all        # то же самое (синоним)
#   ./start.sh build      # только инкрементальная сборка
#   ./start.sh run        # только запуск run_afterburner.C (без сборки)
#   ./start.sh clean      # чистая сборка (удаляет build/)
#   ./start.sh rebuild    # make clean && make
#
# Параметры запуска (inputFile, dmin, p0, metricMode, outputBase)
# задаются прямо в run_afterburner.C.
set -e

ROOT_DIR=/home/oem/WORKLAST/AAMCC-newDataHandle
GEANT4_INSTALL=$ROOT_DIR/geant4-10.4-install
MCINI=$ROOT_DIR/mcini
CXX=/home/oem/miniforge3/bin/g++
BUILD_DIR=$ROOT_DIR/build
JOBS=$(nproc)

usage() {
  echo "Использование: $0 [all|build|run|clean|rebuild]"
  echo "  all      (по умолчанию) — собрать и запустить"
  echo "  build    — инкрементальная сборка"
  echo "  run      — только запуск run_afterburner.C через ROOT"
  echo "  clean    — чистая сборка (удаляет build/)"
  echo "  rebuild  — make clean && make"
  exit 1
}

source "$GEANT4_INSTALL/bin/geant4.sh"
export MCINI="$MCINI"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$GEANT4_INSTALL/lib:$MCINI/build"

do_build() {
  cd "$BUILD_DIR"
  make -j"$JOBS"
  echo "Сборка готова: $BUILD_DIR/GRATE"
}

do_clean() {
  rm -rf "$BUILD_DIR"
  mkdir -p "$BUILD_DIR"
  cd "$BUILD_DIR"
  cmake -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_STANDARD=17 \
        -DCMAKE_CXX_COMPILER="$CXX" \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DCMAKE_MAKE_PROGRAM=/usr/bin/make \
        -DGeant4_DIR="$GEANT4_INSTALL/lib/cmake/Geant4" \
        ..
  make -j"$JOBS"
  echo "Чистая сборка готова: $BUILD_DIR/GRATE"
}

do_rebuild() {
  cd "$BUILD_DIR"
  make clean
  make -j"$JOBS"
  echo "Пересборка готова: $BUILD_DIR/GRATE"
}

do_run() {
  cd "$ROOT_DIR"
  root -l -b -q run_afterburner.C
}

case "${1:-all}" in
  all)
    do_build
    do_run
    ;;
  build)
    do_build
    ;;
  run)
    do_run
    ;;
  clean)
    do_clean
    ;;
  rebuild)
    do_rebuild
    ;;
  *)
    usage
    ;;
esac
