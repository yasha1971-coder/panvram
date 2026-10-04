"""Windows cut straight from FASTA bytes (independent of panvram's own tables): the comparison side of the tests."""
COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")
TOK = bytes(0 if c in b"Aa" else 1 if c in b"Cc" else 2 if c in b"Gg" else 3 if c in b"Tt" else 4 for c in range(256))


class FastaIndex:
    """name -> (offset of the first base in the file, line width, length) by scanning the headers."""

    def __init__(self, fa):
        self.fa, self.rec = fa, {}
        h = 0 if fa[:1] == b">" else -1
        while h >= 0:
            e = fa.find(b"\n", h)
            hdr = fa[h + 1:e]
            name = hdr.split(b" ")[0].split(b"\t")[0].decode()
            nxt = fa.find(b"\n>", e)
            end = len(fa) if nxt < 0 else nxt + 1
            region = end - (e + 1)
            lw = fa.find(b"\n", e + 1) - (e + 1) if region else 0
            n = 0
            if region:
                k, rem = divmod(region, lw + 1)
                n = k * lw + (rem - 1 if rem else 0)
            self.rec[name] = (e + 1, lw, n)
            h = -1 if nxt < 0 else nxt + 1

    def cut(self, name, s, W):
        off, lw, n = self.rec[name]
        assert 0 <= s and s + W <= n, (name, s, W, n)
        if W == 0:
            return b""
        a = off + s + s // lw
        b = off + (s + W - 1) + (s + W - 1) // lw + 1
        out = self.fa[a:b].replace(b"\n", b"")
        assert len(out) == W
        return out


def transform(b, rc=False, tokens=False):
    if rc:
        b = b.translate(COMP)[::-1]
    if tokens:
        b = b.translate(TOK)
    return b
