"""Interrupted numerical campaigns reuse only complete source-bound parts."""
import json

import pytest

import gaussian_dynamics.v307_benchmark as benchmark


def test_case_checkpoint_reuses_exact_source_and_runtime_only(tmp_path,monkeypatch):
    calls=[]
    def compute(arguments):
        calls.append(arguments)
        return {"measured":len(calls)}
    monkeypatch.setattr(benchmark,"_frame_v307",compute)
    task=("frame",("fixture",.1,.3,False,False),"a"*40,str(tmp_path))
    assert benchmark._campaign_part_v307(task) == {"measured":1}
    assert benchmark._campaign_part_v307(task) == {"measured":1}
    assert len(calls) == 1
    changed=(*task[:2],"b"*40,task[3])
    assert benchmark._campaign_part_v307(changed) == {"measured":2}
    monkeypatch.setattr(benchmark.platform,"python_version",lambda:"different-runtime")
    assert benchmark._campaign_part_v307(task) == {"measured":3}


def test_corrupt_or_partial_checkpoint_fails_closed(tmp_path,monkeypatch):
    monkeypatch.setattr(benchmark,"_frame_v307",lambda args:{"measured":1})
    task=("frame",("fixture",.1,.3,False,False),"a"*40,str(tmp_path))
    benchmark._campaign_part_v307(task)
    path=next(tmp_path.glob('*.json'))
    data=json.loads(path.read_text());data['result']['measured']=2;path.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='integrity or provenance'):
        benchmark._campaign_part_v307(task)


def test_failed_case_never_publishes_a_cache_entry(tmp_path,monkeypatch):
    def fail(args):
        raise ValueError('failed solve')
    monkeypatch.setattr(benchmark,"_frame_v307",fail)
    with pytest.raises(ValueError,match='failed solve'):
        benchmark._campaign_part_v307(("frame",(),"a"*40,str(tmp_path)))
    assert not list(tmp_path.iterdir())
