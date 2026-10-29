"""Read SOURMASH-MANIFEST.csv of a remote sourmash zip via HTTP range requests; print summary + sample rows."""
import sys, io, zipfile, csv, collections, requests

class HTTPFile(io.RawIOBase):
    def __init__(self, url):
        self.url = url; self.pos = 0
        r = requests.head(url, allow_redirects=True, timeout=30); self.size = int(r.headers["content-length"])
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]; return self.pos
    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0: return 0
        r = requests.get(self.url, headers={"Range": f"bytes={self.pos}-{self.pos+n-1}"}, timeout=120)
        d = r.content; b[:len(d)] = d; self.pos += len(d); return len(d)

for url in sys.argv[1:]:
    z = zipfile.ZipFile(io.BufferedReader(HTTPFile(url), buffer_size=1 << 20))
    mf = z.read("SOURMASH-MANIFEST.csv").decode()
    lines = mf.splitlines(); rows = list(csv.DictReader(lines[1:]))
    c = collections.Counter((r["ksize"], r["moltype"], r["scaled"], r["num"], r["with_abundance"]) for r in rows)
    nh = sorted(int(r["n_hashes"]) for r in rows)
    print(f"== {url.split('/')[-1]}: {len(rows)} sketches; (ksize,moltype,scaled,num,abund) {dict(c)}")
    print(f"   n_hashes min/p10/median/max: {nh[0]}/{nh[len(nh)//10]}/{nh[len(nh)//2]}/{nh[-1]}; zero: {sum(1 for x in nh if x==0)}; <10: {sum(1 for x in nh if x<10)}")
    for r in rows[:3]: print("   name:", r["name"][:110])
