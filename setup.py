"""Build a platform shared library for ctypes; no Python extension entry point."""

import os
from pathlib import Path
import shutil
import sys

from setuptools import Extension, setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.build_ext import build_ext
from setuptools.errors import PlatformError


class BuildNativeLibrary(build_ext):
    def run(self):
        if sys.platform not in ("linux", "win32", "darwin"):
            raise PlatformError("smoothing_on_sphere supports Linux, Windows and macOS.")
        super().run()

    def get_ext_filename(self, fullname):
        return {"linux": "libsmoothing_on_sphere.so", "win32": "smoothing_on_sphere.dll",
                "darwin": "libsmoothing_on_sphere.dylib"}[sys.platform]

    def get_export_symbols(self, ext):
        return []

    def get_libraries(self, ext):
        return ext.libraries

    def build_extensions(self):
        prefix = os.environ.get("LIBOMP_PREFIX")
        if self.compiler.compiler_type == "msvc":
            if not self.compiler.initialized:
                self.compiler.initialize()
            # Clang accepts the existing unsigned OpenMP loop indices.
            clang = shutil.which("clang-cl")
            if not clang or not prefix:
                raise PlatformError("Windows source builds require clang-cl, MSVC Build Tools, "
                                    "and LIBOMP_PREFIX pointing to an LLVM OpenMP prefix "
                                    "containing lib/libomp.lib and bin/libomp.dll.")
            self.compiler.cc = clang
            compile_args = ["/O2", "/std:c++14", "/EHsc", "/W4", "/D_USE_MATH_DEFINES",
                            "/clang:-fopenmp=libomp"]
            link_args = ["/LIBPATH:" + str(Path(prefix) / "lib"), "libomp.lib"]
        else:
            compile_args = ["-O3", "-std=c++11", "-Wall", "-Wextra", "-Wpedantic",
                            "-fvisibility=hidden"]
            if sys.platform == "darwin":
                if not prefix:
                    # Native Homebrew installations; explicit prefix takes precedence.
                    candidates = (Path("/opt/homebrew/opt/libomp"), Path("/usr/local/opt/libomp"))
                    prefix = next((str(p) for p in candidates if (p / "include/omp.h").is_file()), None)
                if not prefix:
                    raise PlatformError("macOS source builds require libomp; set LIBOMP_PREFIX "
                                        "to its installation prefix.")
                compile_args += ["-Xpreprocessor", "-fopenmp", "-I" + str(Path(prefix) / "include")]
                link_args = ["-L" + str(Path(prefix) / "lib"), "-lomp",
                             "-Wl,-rpath," + str(Path(prefix) / "lib")]
                for name in ("linker_so", "linker_so_cxx"):
                    command = getattr(self.compiler, name, None)
                    if command:
                        setattr(self.compiler, name, ["-dynamiclib" if part == "-bundle" else part
                                for part in command if part not in ("-undefined", "dynamic_lookup")])
            else:
                # Release wheels explicitly select LLVM. Local Linux source
                # builds keep the compiler's default runtime (including GCC).
                runtime = os.environ.get("SMOOTHING_OPENMP_RUNTIME", "default")
                if runtime not in ("default", "llvm"):
                    raise PlatformError("SMOOTHING_OPENMP_RUNTIME must be default or llvm.")
                openmp_flag = "-fopenmp=libomp" if runtime == "llvm" else "-fopenmp"
                compile_args += [openmp_flag]
                link_args = [openmp_flag, "-Wl,--no-undefined"]
                if prefix:
                    compile_args += ["-I" + str(Path(prefix) / "include")]
                    link_args += ["-L" + str(Path(prefix) / "lib"),
                                  "-Wl,-rpath," + str(Path(prefix) / "lib")]
        for extension in self.extensions:
            extension.extra_compile_args = compile_args
            extension.extra_link_args = link_args
        super().build_extensions()


class BuildPlatformWheel(bdist_wheel):
    def get_tag(self):
        _, _, platform = super().get_tag()
        return f"cp{sys.version_info.major}{sys.version_info.minor}", "none", platform


setup(
    ext_modules=[Extension(
        "libsmoothing_on_sphere",
        sources=["src/cpp/CC_smoothing_on_sphere_python_lib.cc"],
        include_dirs=["src/cpp"],
        depends=["src/cpp/CU_smoothing_on_sphere_code.cc", "src/cpp/CU_kdtree_with_index.cc"],
        language="c++",
    )],
    cmdclass={"build_ext": BuildNativeLibrary, "bdist_wheel": BuildPlatformWheel},
)
