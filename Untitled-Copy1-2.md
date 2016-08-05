

```python
##simple alignment
def simp(p,t):
    count=0
    align=0
    n=0
    pos=[]
    for i in range(len(t)-len(p)+1):
        match=True
        n+=1
        for j in range(len(p)):
            count+=1
            if p[j]!=t[i+j]:
                match=False
                align+=1
                break
        if match:
            pos.append(i)
            align1=align+1
    return pos,count,align1,n
    
```


```python
p = 'CTGT'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as

```


```python
t

```




    'AAAAAAAAAACTGTAAAAAAAAAACTTTAAAAAAAAAACGGGAAAAAAAAAA'




```python
simp(ten_as,t)
```




    ([0, 14, 28, 42], 214, 40, 43)




```python
p1="GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
p2="GGCGCGGTGGCTCACGCCTGTAAT"
```


```python
seq=""
fh=open("/Users/daxinggao/Downloads/chr1.GRCh38.excerpt.fasta")
for line in fh:
    if not line.startswith(">"):
        seq+=line.rstrip()

```


```python
simp(p1,seq)

```




    ([56922], 984143, 56923, 799954)




```python
def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    count=0
    align=0
    n=0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        count+=1
        mismatched = False
        n+=1
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                align+=1
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
            align1=align+1
        i += shift
    return occurrences,count,align1,n
```


```python

import unittest


def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]



```


```python
p1_bm=BoyerMoore(p1)
boyer_moore(p1,p1_bm,seq)

```




    ([56922], 127974, 9381, 127974)




```python
import bisect


class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

```


```python
t = 'to-morrow and to-morrow and to-mdrrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
ind=Index(t,10)
ind.query(p)
```




    [0, 14]




```python
def appro(p,t,n):
    matches=[]
    seg=3
    ind=Index(t,8)
    hit=0
    for i in range(seg):
        start=i*8
        end=i*8+8
        cand=ind.query(p[start:end])
        hit+=len(cand)
        ##print(start,end)
        if(len(cand)>0):
            ##print(cand)
            for pos in cand:
                mismatch=0 ## count mismatch
                if pos+len(p)-start<len(t):
                    for j in range(-start,len(p)-start,1):
                        if p[start+j]!=t[pos+j]:
                            mismatch+=1
                if mismatch<=n:
                        matches.append(pos-start)
    return set(matches),len(matches),len(set(matches)),hit


                    
        
    
    
```


```python
appro(p1,seq,2)
```




    ({56922,
      147558,
      160162,
      160729,
      191452,
      364263,
      429299,
      465647,
      657496,
      717706,
      724927},
     25,
     11)




```python
import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
```


```python
def appro2(p,t,n):
    ind2=SubseqIndex(t,8,3)
    matches=[]
    hit=0
    for i in range(3):
        cand2=ind2.query(p[i:])
        hit+=len(cand2)
        for pos in cand2:
            mismatch=0
            for j in range(-i,len(p)-i,1):
                if p[j+i]!=t[pos+j]:
                        mismatch+=1
            if mismatch<=n:matches.append(pos-i)
    return matches,len(set(matches)),len(matches),hit

```


```python
appro2(p2,seq,2)
```




    ([56922,
      84641,
      147558,
      191452,
      262042,
      273669,
      364263,
      465647,
      635931,
      657496,
      681737,
      717706,
      747359,
      56922,
      84641,
      147558,
      160162,
      160729,
      262042,
      273669,
      364263,
      421221,
      429299,
      465647,
      551134,
      635931,
      657496,
      681737,
      717706,
      724927,
      56922,
      160729,
      191452,
      262042,
      364263,
      429299,
      657496,
      717706,
      724927,
      747359],
     19,
     40,
     79)




```python
appro2(p,t,2)
```




    ([0, 0, 14, 0, 14], 2, 5, 5)




```python
def approximate_match_subseq(p, t, n, ival):
    segment_length = int(round(len(p) / (n+1)))
    all_matches = set()
    p_idx = SubseqIndex(t, segment_length, ival)
    idx_hits = 0
    for i in range(n+1):
        start = i
        matches = p_idx.query(p[start:])
        
        # Extend matching segments to see if whole p matches
        for m in matches:
            idx_hits += 1
            if m < start or m-start+len(p) > len(t):
                continue
            
            mismatches = 0
            
            for j in range(0, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches), idx_hits

```


```python
approximate_match_subseq(p2,seq,2,3)
```




    ([84641,
      160162,
      273669,
      147558,
      364263,
      421221,
      681737,
      717706,
      724927,
      465647,
      429299,
      657496,
      160729,
      56922,
      635931,
      191452,
      262042,
      551134,
      747359],
     79)




```python
appro(p2,seq,2)

```




    ({56922,
      84641,
      147558,
      160162,
      160729,
      191452,
      262042,
      273669,
      364263,
      421221,
      429299,
      465647,
      551134,
      635931,
      657496,
      681737,
      717706,
      724927,
      747359},
     40,
     19,
     90)


