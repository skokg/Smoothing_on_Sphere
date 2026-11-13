"""Numerical/API checks for an installed package; no graphical dependencies.

Run explicitly with: python tests/test_smoothing_on_sphere.py
Disk caches use temporary directories. This file has not been executed during
packaging preparation; builds and tests are performed by the user or CI.
"""
import copy
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
import smoothing_on_sphere as smoothing


def reference_smoothing(lat, lon, area, field, radius):
    """Independent all-pairs spherical-cap average, away from cutoff ties."""
    area = np.asarray(area, dtype=np.float64)
    field = np.asarray(field, dtype=np.float64)
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    hav = np.sin(dlat / 2) ** 2 + (
        np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2)
    distances = 2 * 6371000.0 * np.arcsin(np.sqrt(np.clip(hav, 0, 1)))
    result = np.zeros(field.shape, dtype=np.float64)
    for i in range(field.size):
        included = distances[i] <= radius
        weight = area[included].sum()
        if area[i] > 0 and weight > 0:
            result[i] = np.dot(area[included], field[included]) / weight
    return result


def cache_variants():
    return (
        (smoothing.generate_kdtree_cache_data_in_memory,
         smoothing.save_kdtree_cache,
         smoothing.load_kdtree_cache,
         smoothing.smooth_field_using_kdtree_cache_data,
         smoothing.smooth_multiple_fields_simultaneously_using_kdtree_cache_data,
         smoothing.free_kdtree_cache_data),
        (smoothing.generate_overlap_cache_data_in_memory,
         smoothing.save_overlap_cache,
         smoothing.load_overlap_cache,
         smoothing.smooth_field_using_overlap_detection,
         smoothing.smooth_multiple_fields_simultaneously_using_overlap_detection,
         smoothing.free_overlap_cache_data))

class SmoothingTests(unittest.TestCase):
    # These tests intentionally use small inputs: even the larger reference
    # problem has only 320 points, and every reference is an all-pairs calculation.
    def setUp(self):
        self.lat = np.zeros(4)
        self.lon = np.array([0., 1., 3., 6.])
        self.area = np.array([1., 3., 2., 4.])
        self.field = np.array([2., 10., 20., 30.])
        self.radius = 150000.
        self.tree = smoothing.construct_KdTree(self.lat, self.lon, 1)

    def tearDown(self):
        if self.tree is not None:
            smoothing.free_KdTree_memory(self.tree)

    def direct(self, field=None, area=None, threads=1, radius=None):
        return smoothing.smooth_field_using_KdTree(
            self.radius if radius is None else radius, self.tree, self.lat, self.lon,
            self.area if area is None else area, self.field if field is None else field,
            threads, False)

    def test_native_abi(self):
        self.assertEqual(smoothing.libc.smoothing_wrapper_abi_version(), 3)

    def test_known_area_weighted_average(self):
        actual = self.direct()
        np.testing.assert_allclose(actual, [8., 8., 20., 30.], rtol=1e-12, atol=1e-12)
        self.assertEqual(actual.dtype, np.dtype('float64'))

    def test_constant_field(self):
        np.testing.assert_allclose(self.direct(np.full(4, 7.)), 7.)

    def test_small_radius_preserves_isolated_values(self):
        np.testing.assert_allclose(self.direct(radius=100.), self.field)

    def test_missing_cells(self):
        area = np.array([1., 0., 2., 4.])
        field = self.field.copy()
        field[1] = 1e20  # Finite placeholder; zero area excludes its contribution.
        np.testing.assert_allclose(self.direct(field, area), [2., 0., 20., 30.])

    def test_spherical_reference_dateline_and_poles(self):
        lat = np.array([0., 0., 35., 35., 89., 89., -40., -40.])
        lon = np.array([179.5, -179.5, 20., 21., 0., 90., -80., -81.])
        area = np.arange(1., 9.)
        field = np.array([1., 8., 2., 7., 3., 6., 4., 5.])
        tree = smoothing.construct_KdTree(lat, lon, 2)
        try:
            actual = smoothing.smooth_field_using_KdTree(
                200000., tree, lat, lon, area, field, 2, False)
            np.testing.assert_allclose(actual, reference_smoothing(
                lat, lon, area, field, 200000.), rtol=1e-12, atol=1e-12)
        finally:
            smoothing.free_KdTree_memory(tree)

    def test_direct_multiple_fields_and_threads(self):
        fields = np.array([self.field, 2 * self.field + 1])
        expected = np.array([reference_smoothing(
            self.lat, self.lon, self.area, f, self.radius) for f in fields])
        for threads in (1, 2):
            actual = smoothing.smooth_multiple_fields_simultaneously_using_KdTree(
                self.radius, self.tree, self.lat, self.lon, self.area, fields, threads, False)
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(self.direct(threads=threads), expected[0])

    def check_cache(self, cache, single, multiple, radius):
        for area in (self.area, np.array([1., 0., 2., 4.])):
            fields = np.array([self.field, 2 * self.field + 1])
            expected = np.array([reference_smoothing(
                self.lat, self.lon, area, f, radius) for f in fields])
            for threads in (1, 2):
                np.testing.assert_allclose(single(area, fields[0], cache, threads, False),
                                           expected[0], rtol=1e-12, atol=1e-12)
                np.testing.assert_allclose(multiple(area, fields, cache, threads, False),
                                           expected, rtol=1e-12, atol=1e-12)

    def test_in_memory_kdtree_cache(self):
        cache = smoothing.generate_kdtree_cache_data_in_memory(
            self.lat, self.lon, self.radius, 2, False)
        try:
            self.check_cache(cache, smoothing.smooth_field_using_kdtree_cache_data,
                             smoothing.smooth_multiple_fields_simultaneously_using_kdtree_cache_data,
                             self.radius)
        finally:
            smoothing.free_kdtree_cache_data(cache)
        with self.assertRaises(ValueError):
            smoothing.free_kdtree_cache_data(cache)

    def test_disk_round_trip_for_each_radius(self):
        for generate, write, load, single, multiple, free in cache_variants():
            for radius in (self.radius, 400000.):
                with self.subTest(method=generate.__name__, radius=radius), tempfile.TemporaryDirectory() as folder:
                    cache = generate(self.lat, self.lon, radius, 2, False)
                    try:
                        self.check_cache(cache, single, multiple, radius)
                        write(cache, Path(folder) / "custom cache.bin", False)
                    finally:
                        free(cache)
                    files = list(Path(folder).iterdir())
                    self.assertEqual(files, [Path(folder) / "custom cache.bin"])
                    original = files[0].read_bytes()
                    loaded = load(Path(folder) / "custom cache.bin", False)
                    try:
                        self.check_cache(loaded, single, multiple, radius)
                        # A loaded cache can be saved again without regeneration.
                        write(loaded, Path(folder) / "custom cache.bin", False)
                        self.assertEqual(files[0].read_bytes(), original)
                    finally:
                        free(loaded)
                    with self.assertRaises(ValueError):
                        write(loaded, Path(folder) / "custom cache.bin", False)
                    with self.assertRaises(ValueError):
                        free(loaded)

    def test_invalid_inputs(self):
        for radius in (-1., float('nan'), float('inf')):
            with self.subTest(radius=radius), self.assertRaises(ValueError):
                self.direct(radius=radius)
        for threads in (0, -1):
            with self.subTest(threads=threads), self.assertRaises(ValueError):
                self.direct(threads=threads)
        with self.assertRaises(TypeError):
            self.direct(threads=True)
        with self.assertRaises(ValueError):
            self.direct(field=np.array([1., 2.]))
        with self.assertRaises(ValueError):
            self.direct(field=np.array([1., 2., float('nan'), 4.]))
        with self.assertRaises(TypeError):
            self.direct(field=np.ma.array(self.field, mask=False))

    def test_handle_ownership_and_cleanup(self):
        with self.assertRaises(TypeError):
            copy.copy(self.tree)
        with self.assertRaises(TypeError):
            copy.deepcopy(self.tree)
        tree = self.tree
        smoothing.free_KdTree_memory(tree)
        self.tree = None
        with self.assertRaises(ValueError):
            smoothing.free_KdTree_memory(tree)
        with self.assertRaises(ValueError):
            smoothing.smooth_field_using_KdTree(
                self.radius, tree, self.lat, self.lon, self.area, self.field, 1, False)

    def test_missing_cache_reports_native_error(self):
        with tempfile.TemporaryDirectory() as folder:
            for load in (smoothing.load_kdtree_cache,
                         smoothing.load_overlap_cache):
                with self.subTest(loader=load.__name__), self.assertRaises(RuntimeError):
                    load(Path(folder) / "custom cache.bin", False)

    def check_all_methods(self, lat, lon, area, fields, radii, threads=(1, 2)):
        """Compare every smoothing route with an independent numerical oracle."""
        expected = {
            radius: np.array([reference_smoothing(lat, lon, area, field, radius)
                              for field in fields]) for radius in radii}

        def check(single, multiple, cache, radius):
            for count in threads:
                with self.subTest(method=single.__name__, radius=radius, threads=count):
                    actual = single(area, fields[0], cache, count, False)
                    batch = multiple(area, fields, cache, count, False)
                    self.assertEqual(actual.shape, area.shape)
                    self.assertEqual(batch.shape, fields.shape)
                    self.assertEqual(actual.dtype, np.dtype('float64'))
                    self.assertEqual(batch.dtype, np.dtype('float64'))
                    np.testing.assert_allclose(actual, expected[radius][0],
                                               rtol=1e-12, atol=1e-12)
                    np.testing.assert_allclose(batch, expected[radius],
                                               rtol=1e-12, atol=1e-12)

        tree = smoothing.construct_KdTree(lat, lon, 1)
        try:
            for radius in radii:
                for count in threads:
                    with self.subTest(method='direct', radius=radius, threads=count):
                        actual = smoothing.smooth_field_using_KdTree(
                            radius, tree, lat, lon, area, fields[0], count, False)
                        batch = smoothing.smooth_multiple_fields_simultaneously_using_KdTree(
                            radius, tree, lat, lon, area, fields, count, False)
                        np.testing.assert_allclose(actual, expected[radius][0],
                                                   rtol=1e-12, atol=1e-12)
                        np.testing.assert_allclose(batch, expected[radius],
                                                   rtol=1e-12, atol=1e-12)
        finally:
            smoothing.free_KdTree_memory(tree)

        for generate, write, load, single, multiple, free in cache_variants():
            for radius in radii:
                # Separate files per radius also allow testing radii that cap to
                # the same whole-sphere value without filename collisions.
                with tempfile.TemporaryDirectory() as folder:
                    cache = generate(lat, lon, radius, 1, False)
                    try:
                        check(single, multiple, cache, radius)
                        write(cache, Path(folder) / "custom cache.bin", False)
                    finally:
                        free(cache)
                    cache = load(Path(folder) / "custom cache.bin", False)
                    try:
                        check(single, multiple, cache, radius)
                    finally:
                        free(cache)

    def test_larger_reference_across_all_methods(self):
        # Shuffling removes the relationship between input order and geography.
        # More than 256 points also exercises a wider range of cache indices.
        rng = np.random.default_rng(20261112)
        lat, lon = np.meshgrid(np.linspace(-82., 82., 16),
                               np.linspace(-171., 171., 20), indexing='ij')
        order = rng.permutation(lat.size)
        lat, lon = lat.ravel()[order], lon.ravel()[order]
        area = rng.uniform(0.1, 4., lat.size)
        area[::11] = 0.
        fields = rng.normal(size=(3, lat.size))
        self.check_all_methods(lat, lon, area, fields, [2000000., 7000000.])

    def test_dateline_and_poles_across_all_methods(self):
        lat = np.array([0., 0., 89., 89., -89., -89., 90., -90.])
        lon = np.array([179.5, -179.5, 0., 90., -40., 50., 180., -180.])
        area = np.arange(1., 9.)
        fields = np.array([np.arange(8.), np.arange(8.)[::-1]])
        self.check_all_methods(lat, lon, area, fields, [200000.])

    def test_boundary_cases_across_all_methods(self):
        cases = (
            ('one-point', np.array([0.]), np.array([0.]), np.array([2.]),
             np.array([[7.], [-3.]]), [0.]),
            ('duplicates-zero-radius', np.zeros(4), np.array([0., 0., 3., 6.]),
             self.area, np.array([self.field, -self.field]), [0.]),
            ('all-missing', self.lat, self.lon, np.zeros(4),
             np.array([self.field, -self.field]), [self.radius]),
            ('whole-sphere', np.array([0., 0., 90., -90.]),
             np.array([0., 180., 0., 0.]), self.area,
             np.array([self.field, -self.field]),
             [np.pi * 6371000., 2 * np.pi * 6371000.]))
        for name, lat, lon, area, fields, radii in cases:
            with self.subTest(case=name):
                self.check_all_methods(lat, lon, area, fields, radii, threads=(1,))

    def test_array_conversions_across_all_methods(self):
        for dtype in (np.int32, np.float32):
            # Every input is a strided view; fields are also not C contiguous.
            def strided(values):
                values = np.asarray(values, dtype=dtype)
                storage = np.zeros(values.shape[:-1] + (values.shape[-1] * 2,), dtype=dtype)
                storage[..., ::2] = values
                return storage[..., ::2]

            with self.subTest(dtype=dtype):
                inputs = [strided(values) for values in (
                    self.lat, self.lon, self.area, np.array([self.field, -self.field]))]
                originals = [values.copy() for values in inputs]
                self.check_all_methods(*inputs, [self.radius], threads=(1,))
                for actual, original in zip(inputs, originals):
                    np.testing.assert_array_equal(actual, original)

    def test_additional_invalid_arrays_and_counts(self):
        for value, error in (
                (np.array([]), ValueError),
                (np.ones((2, 2)), ValueError),
                (np.array([1., 2., np.inf, 4.]), ValueError),
                (np.array([1., 2., 3., 4.], dtype=complex), TypeError),
                (np.array(['1', '2', '3', '4']), TypeError),
                ([1., 2., 3., 4.], TypeError)):
            with self.subTest(field=repr(value)), self.assertRaises(error):
                self.direct(field=value)
        for threads, error in ((1.5, TypeError), (np.bool_(True), TypeError),
                               (2 ** 31, ValueError)):
            with self.subTest(threads=threads), self.assertRaises(error):
                self.direct(threads=threads)
        with self.assertRaises(ValueError):
            smoothing.construct_KdTree(self.lat, self.lon[:-1], 1)
        with self.assertRaises(ValueError):
            smoothing.smooth_multiple_fields_simultaneously_using_KdTree(
                self.radius, self.tree, self.lat, self.lon, self.area,
                np.ones((2, 3)), 1, False)
        with self.assertRaises(TypeError):
            smoothing.smooth_field_using_kdtree_cache_data(
                self.area, self.field, self.tree, 1, False)

    def test_damaged_disk_caches(self):
        for generate, write, load, single, multiple, free in cache_variants():
            with tempfile.TemporaryDirectory() as folder:
                cache = generate(self.lat, self.lon, self.radius, 1, False)
                try:
                    write(cache, Path(folder) / "custom cache.bin", False)
                finally:
                    free(cache)
                files = list(Path(folder).glob('*.bin'))
                self.assertEqual(len(files), 1)
                path = files[0]
                original = path.read_bytes()
                for label, damaged in (
                        ('empty', b''), ('short-header', original[:4]),
                        ('bad-magic', b'BADMAGIC' + original[8:]),
                        ('truncated-record', original[:-1]),
                        ('trailing-data', original + b'X')):
                    with self.subTest(method=load.__name__, damage=label):
                        path.write_bytes(damaged)
                        cache = None
                        try:
                            with self.assertRaises(RuntimeError):
                                cache = load(Path(folder) / "custom cache.bin", False)
                        finally:
                            if cache is not None:
                                free(cache)
                path.write_bytes(original)
                cache = load(Path(folder) / "custom cache.bin", False)
                try:
                    np.testing.assert_allclose(single(self.area, self.field, cache, 1, False),
                                               [8., 8., 20., 30.])
                finally:
                    free(cache)

    def test_overlap_multiblock_save_and_reload(self):
        # Whole-sphere records include a 300-point first record (more than one
        # 256-index block) followed by records with empty added/removed sets.
        count = 300
        lat = np.zeros(count)
        lon = np.linspace(-180., 180., count, endpoint=False)
        area = 1. + np.arange(count) % 7
        field = np.arange(count, dtype=np.float64)
        radius = np.pi * 6371000.
        expected = np.full(count, np.dot(area, field) / area.sum())
        with tempfile.TemporaryDirectory() as folder:
            cache = smoothing.generate_overlap_cache_data_in_memory(lat, lon, radius, 2, False)
            try:
                np.testing.assert_allclose(smoothing.smooth_field_using_overlap_detection(
                    area, field, cache, 2, False), expected, rtol=1e-12, atol=1e-12)
                smoothing.save_overlap_cache(cache, Path(folder) / "custom cache.bin", False)
            finally:
                smoothing.free_overlap_cache_data(cache)
            path, = Path(folder).iterdir()
            original = path.read_bytes()
            cache = smoothing.load_overlap_cache(Path(folder) / "custom cache.bin", False)
            try:
                np.testing.assert_allclose(smoothing.smooth_field_using_overlap_detection(
                    area, field, cache, 2, False), expected, rtol=1e-12, atol=1e-12)
                smoothing.save_overlap_cache(cache, Path(folder) / "custom cache.bin", False)
                self.assertEqual(path.read_bytes(), original)
            finally:
                smoothing.free_overlap_cache_data(cache)

    def test_cache_load_uses_contents_not_filename(self):
        for generate, write, load, single, multiple, free in cache_variants():
            for radius in (self.radius, 400000.):
                with self.subTest(method=generate.__name__, radius=radius), tempfile.TemporaryDirectory() as folder:
                    filename = Path(folder) / 'chosen filename.dat'
                    cache = generate(self.lat, self.lon, radius, 1, False)
                    try:
                        write(cache, filename=filename, show_output_log=False)
                    finally:
                        free(cache)
                    self.assertEqual(list(Path(folder).iterdir()), [filename])
                    renamed = filename.rename(Path(folder) / 'renamed without extension')
                    cache = load(filename=str(renamed), show_output_log=False)
                    try:
                        self.check_cache(cache, single, multiple, radius)
                    finally:
                        free(cache)

    def test_cache_round_trip_with_unicode_paths(self):
        # Exercise UTF-16 surrogate pairs as well as accented and CJK characters.
        unicode_name = '\u010d\u0161\u017e_\u6e2c\u8a66_\U0001f30d'
        for generate, write, load, single, multiple, free in cache_variants():
            with self.subTest(method=generate.__name__), tempfile.TemporaryDirectory() as folder:
                directory = Path(folder) / ('cache_' + unicode_name)
                directory.mkdir()
                filename = directory / (unicode_name + '.bin')
                cache = generate(self.lat, self.lon, self.radius, 1, False)
                try:
                    write(cache, filename, False)
                finally:
                    free(cache)
                self.assertEqual(list(directory.iterdir()), [filename])
                previous = os.getcwd()
                try:
                    os.chdir(directory)
                    for path in (filename, Path(filename.name)):
                        with self.subTest(path=str(path)):
                            loaded = load(path, False)
                            try:
                                self.check_cache(loaded, single, multiple, self.radius)
                            finally:
                                free(loaded)
                finally:
                    os.chdir(previous)

    def test_single_radius_and_no_construction_files(self):
        for generate, write, load, single, multiple, free in cache_variants():
            for radius in ([self.radius], [self.radius, 400000.], np.array([self.radius])):
                with self.subTest(method=generate.__name__, radius=radius), self.assertRaises(ValueError):
                    generate(self.lat, self.lon, radius, 1, False)
            with tempfile.TemporaryDirectory() as folder:
                previous = os.getcwd()
                try:
                    os.chdir(folder)
                    cache = generate(self.lat, self.lon, self.radius, 1, False)
                    try:
                        self.assertEqual(list(Path('.').iterdir()), [])
                        np.testing.assert_allclose(single(self.area, self.field, cache, 1, False),
                                                   [8., 8., 20., 30.])
                    finally:
                        free(cache)
                finally:
                    os.chdir(previous)

    def test_failed_save_preserves_cache_and_existing_file(self):
        for generate, write, load, single, multiple, free in cache_variants():
            with tempfile.TemporaryDirectory() as folder:
                cache = generate(self.lat, self.lon, self.radius, 1, False)
                try:
                    write(cache, Path(folder) / "custom cache.bin", False)
                    path, = Path(folder).iterdir()
                    original = path.read_bytes()
                    # Portable failure injection at publication, after native writing.
                    with mock.patch.object(smoothing.os, 'replace', side_effect=OSError('simulated write failure')):
                        with self.assertRaises(OSError):
                            write(cache, Path(folder) / "custom cache.bin", False)
                    self.assertEqual(path.read_bytes(), original)
                    self.assertEqual(list(Path(folder).iterdir()), [path])
                    np.testing.assert_allclose(single(self.area, self.field, cache, 1, False),
                                               [8., 8., 20., 30.])
                    with self.assertRaises(TypeError):
                        write(self.tree, Path(folder) / "custom cache.bin", False)
                finally:
                    free(cache)



if __name__ == '__main__':
    unittest.main()
