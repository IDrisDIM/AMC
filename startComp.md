# Инструкция по компиляции AAMCC

## 1. Требования

- **Geant4 10.4** — установлен в `/home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install/`
- **ROOT 6+** — установлен в `/home/oem/miniforge3/`
- **MCini** — лежит в `/home/oem/WORKLAST/AAMCC-newDataHandle/mcini/`
- **CMake 3.0+**, **GCC 7+** (для C++17)

## 2. Активация окружения (каждый новый терминал)

```bash
source /home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install/bin/geant4.sh
export MCINI=/home/oem/WORKLAST/AAMCC-newDataHandle/mcini
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCINI/build
```

Проверка:

```bash
which geant4-config
geant4-config --version    # должно быть 10.4.3
which root
root-config --version      # 6.x
```

## 3. Сборка MCini (один раз, если ещё не собрана)

```bash
cd /home/oem/WORKLAST/AAMCC-newDataHandle/mcini
mkdir -p build && cd build
cmake -DCMAKE_CXX_STANDARD=17 ..
make -j4
```

Проверка:

```bash
ls $MCINI/build/libMcIniData.*
```

## 4. Сборка AAMCC

### 4.1. Чистая сборка (первый раз или после серьёзных правок)

```bash
cd /home/oem/WORKLAST/AAMCC-newDataHandle
rm -rf build
mkdir build && cd build

cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
      -DGeant4_DIR=/home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install/lib/cmake/Geant4 \
      ..

make -j4
```

### 4.2. Пересборка после мелких правок

Если меняли только `.cc` или `.hh`:

```bash
cd /home/oem/WORKLAST/AAMCC-newDataHandle/build
make -j4
```

CMake сам определит, какие файлы пересобрать. Обычно занимает 1–3 минуты.

### 4.3. Полная пересборка (если что-то пошло не так)

```bash
cd /home/oem/WORKLAST/AAMCC-newDataHandle/build
make clean
make -j4
```

## 5. Проверка результата

```bash
ls -la /home/oem/WORKLAST/AAMCC-newDataHandle/build/GRATE
```

Должен быть исполняемый файл. Проверка запуска:

```bash
cd /home/oem/WORKLAST/AAMCC-newDataHandle/build
./GRATE
```

Программа начнёт задавать вопросы. Прервите `Ctrl+C`.

## 6. Типичные ошибки и что делать

### 6.1. `Could not find Geant4Config.cmake`

Не активировано окружение Geant4. Выполните `source .../geant4.sh` и пересоберите.

### 6.2. `fatal error: URun.h: No such file or directory`

Не подключён MCini. Проверьте, что `$MCINI` установлен и `libMcIniData.so` собран в `$MCINI/build`.

### 6.3. `'EnableImplicitMT' is not a member of 'ROOT'`

Закомментируйте `ROOT::EnableImplicitMT()` в `src/AAMCCWriter.cc`.

### 6.4. `error: 'make_unique' is not a member of 'std'`

Не подхватился C++17. Проверьте, что в `CMakeLists.txt` есть:

```cmake
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
```

И в `TGlauber/CMakeLists.txt` тоже.

### 6.5. `undefined reference to ...`

Не хватает библиотеки при линковке. Проверьте, что в `CMakeLists.txt` есть `target_link_libraries(GRATE ... TGlauberLib McIniData)`.

### 6.6. Сборка падает с `-Werror`

Уберите `-Werror` из флагов или добавьте `-Wno-error`.

## 7. Флаги, которые уже проверены и работают

```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
      -DGeant4_DIR=/home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install/lib/cmake/Geant4 \
      ..
```

- `Release` — включает `-O3 -DNDEBUG`, быстрее в 2–3 раза.
- `CMAKE_CXX_COMPILER=/usr/bin/g++` — использует системный компилятор, а не из miniforge3 (иначе возможны конфликты заголовков).
- `Geant4_DIR` — явно указывает путь к конфигурации Geant4 10.4.

## 8. Скрипт для быстрого запуска

Сохраните как `rebuild.sh` в корне AAMCC:

```bash
#!/bin/bash
set -e

GEANT4_INSTALL=/home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install
ROOT_DIR=/home/oem/WORKLAST/AAMCC-newDataHandle

source $GEANT4_INSTALL/bin/geant4.sh
export MCINI=$ROOT_DIR/mcini
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GEANT4_INSTALL/lib:$MCINI/build

cd $ROOT_DIR/build

if [ "$1" == "clean" ]; then
    rm -rf *
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_STANDARD=17 \
          -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
          -DGeant4_DIR=$GEANT4_INSTALL/lib/cmake/Geant4 \
          ..
fi

make -j4
echo "Готово: $ROOT_DIR/build/GRATE"
```

Использование:

```bash
chmod +x rebuild.sh

./rebuild.sh           # обычная пересборка
./rebuild.sh clean     # полная чистая сборка
```

## 9. Полный цикл с нуля (если всё сломалось)

```bash
# 1. Активировать окружение
source /home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install/bin/geant4.sh
export MCINI=/home/oem/WORKLAST/AAMCC-newDataHandle/mcini
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install/lib:$MCINI/build

# 2. Пересобрать MCini (если нужно)
cd $MCINI && rm -rf build && mkdir build && cd build
cmake -DCMAKE_CXX_STANDARD=17 .. && make -j4

# 3. Пересобрать AAMCC
cd /home/oem/WORKLAST/AAMCC-newDataHandle/build
rm -rf *
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=17 \
      -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
      -DGeant4_DIR=/home/oem/WORKLAST/AAMCC-newDataHandle/geant4-10.4-install/lib/cmake/Geant4 \
      ..
make -j4

# 4. Проверить
ls -la GRATE
```

## 10. Время сборки

| Что | Время |
|---|---|
| MCini | 1–2 мин |
| AAMCC с нуля | 5–10 мин |
| AAMCC после мелких правок | 1–3 мин |
| `make clean && make` | 5–10 мин |