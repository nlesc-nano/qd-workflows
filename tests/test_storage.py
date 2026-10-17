"""Storage of a cluster account (qdw campaign status): quota tools of NHR and Hyperion."""
from qd_workflows import storage

NHR = """filesystem,label,allocation,blockUsage,blockSoft,blockHard,blockGrace,inodeUsage,inodeSoft,inodeHard,inodeGrace,lastUpdated
Home,/home/becivinf,users,32471203840,42949672960,107374182400,none,228604,0,0,none,2026-10-17T23:00:13
Work,/scratch/projects/bep00157,projects,2350333722624,13194139533312,131941395333120,none,95376,1000000,6000000,none,2026-10-17T23:00:13
Work,/scratch/usr/becivinf,users,2362845986816,3298534883328,32985348833280,none,237001,250000,1500000,none,2026-10-17T23:00:13
"""
HYP_HOME = """      user/group     ||           size          ||    chunk files
     name     |  id  ||    used    |    hard    ||  used   |  hard
--------------|------||------------|------------||---------|---------
      iinfante|  5730||   20.11 GiB|   50.00 GiB||   282572|   300000
"""


def test_nhr_areas_and_warnings():
    res = storage.storage({"storage": [{"tool": "nhr", "command": "q"}]}, run=lambda c: NHR)
    a = {x["path"]: x for x in res["areas"]}
    assert a["/scratch/usr/becivinf"]["files_frac"] == 0.948 and a["/scratch/projects/bep00157"]["bytes_frac"] == 0.178
    keys = {p["key"]: p["level"] for p in res["problems"]}
    assert keys == {"storage_files_warn:/scratch/usr/becivinf": "warn"}      # 94.8 %: warning, not yet 95 %
    assert a["/home/becivinf"]["files_frac"] is None                         # no file limit on home


def test_hyperion_home_files():
    res = storage.storage({"storage": [{"tool": "beegfs", "command": "getquota-home", "name": "/dipc/iinfante"}]},
                          run=lambda c: HYP_HOME)
    (a,) = res["areas"]
    assert a["bytes"] == int(20.11 * 2 ** 30) and a["files"] == 282572 and a["files_frac"] == 0.942
    assert [p["level"] for p in res["problems"]] == ["warn"]


def test_tool_missing():
    res = storage.storage({"storage": [{"tool": "nhr", "command": "q"}]}, run=lambda c: None)
    assert res["areas"] == [] and res["problems"][0]["level"] == "warn"
