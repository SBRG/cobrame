from __future__ import division, absolute_import, print_function

import pytest

from os.path import dirname, join, abspath

from cobrame.io.jsonme import (load_json_me, save_full_me_model_json,
                               load_full_me_model_json)
from ecolime import build_me_model
from ecolime.util.me_model_comparison import find_me_model_difference

current_dir = dirname(abspath(build_me_model.__file__))
models_dir = current_dir + '/me_models/'

del dirname, abspath

test_model = build_me_model.return_me_model()


def test_full_json_dumping():
    save_full_me_model_json(test_model, 'test_json_dump.json')
    model = load_full_me_model_json('test_json_dump.json')
    difference = find_me_model_difference(model, test_model, 1e-6)
    print(difference)
    assert (len(difference) == 0)


def test_ecoli_build():
    benchmark_model = load_json_me(join(models_dir, 'iJL1678b_benchmark.json'))
    difference = find_me_model_difference(benchmark_model, test_model, 1e-6)
    print('-----------------------Difference----------------------')
    print(difference)
    assert (len(difference) == 0)

if __name__ == '__main__':
    test_ecoli_build()
