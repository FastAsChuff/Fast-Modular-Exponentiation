uint64_t tou64mg(uint64_t a, uint64_t n, uint64_t *twoto64modn) {
  // Converts a to Montgomery form mod n with r = 2^64, n odd, and sets twoto64modn if requested.
  // Recommeded usage...
  // uint64_t twoto64modn = 0;
  // uint64_t ar = tou64mg(a, n, &twoto64modn);
  // uint64_t br = tou64mg(b, n, &twoto64modn); ...etc.
  uint64_t twoto64modnlocal;
  if (n == 0) return 0;
  if ((n & 1) == 0) return 0;
  if (twoto64modn == NULL) {
    twoto64modnlocal = (((unsigned __int128)1) << 64) % n; 
  } else {
    if (*twoto64modn != 0) {
      twoto64modnlocal = *twoto64modn;
    } else {
      twoto64modnlocal = (((unsigned __int128)1) << 64) % n; 
      *twoto64modn = twoto64modnlocal;
    }
  }
  return ((unsigned __int128)twoto64modnlocal*a) % n;
}

uint64_t fromu64mg(uint64_t ar, uint64_t n, uint64_t ninv, uint64_t twoto64modn) {
  // Converts ar out of Montgomery form mod n with r = 2^64, n odd, ninv = n^(-1) mod 2^64.
  uint64_t m = -ar*ninv;
  unsigned __int128 mn128 = (unsigned __int128)m*n; 
  unsigned __int128 sum = ar + mn128;
  unsigned __int128 max128 = (ar > mn128 ? ar : mn128);   
  uint64_t sumu64 = sum >> 64;    
  return (sumu64 + (sum < max128)*twoto64modn) % n;       
}

uint64_t modprodu64mg(uint64_t ar, uint64_t br, uint64_t n, uint64_t ninv, uint64_t twoto64modn) {
  // n must be odd.
  // ar and br must be in Montgomery form with r = 2^64.
  // Returned value is congruent to abr mod n. (NOT necessarily equal to abr mod n.)
  // ab + mn = 0 mod 2^64
  // m = (-ab)ninv mod 2^64
  unsigned __int128 prod128 = (unsigned __int128)ar*br; // 0 <= prod128 <= 2^128 - 2^65 + 1 
  uint64_t m = prod128 & 0xffffffffffffffffULL;         
  m *= -ninv;                                           // 0 <= m < 2^64
  unsigned __int128 mn128 = (unsigned __int128)m*n;     // 0 <= mn128 <= n*(2^64 - 1) && mn128 = -arbr mod 2^64
  unsigned __int128 sum = prod128 + mn128;              // 0 <= sum <= 2^128 + (n - 3)2^64
  unsigned __int128 max128 = (prod128 > mn128 ? prod128 : mn128);
  uint64_t sumu64 = sum >> 64;                          // 0 <= sum >> 64 <= 2^64 + (n - 3)
  // twoto64modn = 2^64 - xn for some x >= 1, so if sum overflowed, 
  // 0 <= sumu64 + twoto64modn <= (n - 3) + 2^64 - xn = 2^64 - (x-1)n - 3 < 2^64
  return sumu64 + (sum < max128)*twoto64modn;
}

uint64_t modpowu64mg(uint64_t ar, uint64_t e, uint64_t n, uint64_t ninv, uint64_t twoto64modn) {
// Returns (a^e)r mod n for odd n, r = 2^64.
  if (n < 2) return 0;
  if ((n & 1) == 0) return 0;
  if (ar < 2) return ar;
  uint64_t res = twoto64modn;
  uint64_t sq = ar;
  while (1) {
    res = ((e & 1ULL) ? modprodu64mg(res, sq, n, ninv, twoto64modn) : res);
    e >>= 1;
    if (e == 0) break;
    sq = modprodu64mg(sq, sq, n, ninv, twoto64modn);
  }
  return res % n;
}

uint32_t modpowu64b(uint32_t a, uint64_t e, uint32_t n) {
  uint32_t res = 1;
  uint32_t sq = a;
  while (1) {
    if (e & 1ULL) res = ((uint64_t)res * sq) % n;
    e >>= 1;
    if (e == 0) break;
    sq = ((uint64_t)sq*sq) % n;
  }
  return res;
}

uint64_t modpowu64general(uint64_t a, uint64_t e, uint64_t n) {
// Returns a^e mod n.
  if (n < 2) return 0;
  if (a < 2) return a;
  if (n <= 0xffffffff) return modpowu64b(a % n, e, n);
  uint64_t res = 1;
  uint64_t sq = a;
  while (1) {
    res = ((e & 1ULL) ? ((unsigned __int128)res * sq) % n : res);
    e >>= 1;
    if (e == 0) break;
    sq = ((unsigned __int128)sq*sq) % n;
  }
  return res;
}

uint64_t modpowu64oddn(uint64_t a, uint64_t e, uint64_t n) {
// Returns (a^e) mod odd n.
  if (n < 2ULL) return 0;
  if (a < 2ULL) return a;
  uint64_t twoto64modn = 0;
  uint64_t ar = tou64mg(a, n, &twoto64modn);
#ifdef USEHACKERSDELIGHTFNS
  uint64_t ninv = modinv64x(n); // Find this here https://github.com/FastAsChuff/Fast-Modular-Inverse-Modulo-Powers-Of-2
#else
  uint64_t ninv = modinv64(n); // Find this here https://github.com/FastAsChuff/Fast-Modular-Inverse-Modulo-Powers-Of-2
#endif
  uint64_t res = modpowu64mg(ar, e, n, ninv, twoto64modn);
  return fromu64mg(res, n, ninv, twoto64modn);
}


uint64_t modpowu64(uint64_t a, uint64_t e, uint64_t n) {
// Returns (a^e) mod n.
  if (n < 2ULL) return 0;
  if (a < 2ULL) return a;
  if (n & 1ULL) {
    return modpowu64oddn(a, e, n);
  } else {
    uint64_t nover2 = n/2ULL;
    if (nover2 & 1ULL) { 
      uint64_t resmodnover2 = modpowu64oddn(a, e, nover2);
      return ((resmodnover2 - a) & 1ULL)*nover2 + resmodnover2;
    } else {
      uint64_t nover4 = n/4ULL;
      if ((nover4 & 1ULL) && (a & 1ULL)) { 
        uint64_t resmodnover4 = modpowu64oddn(a, e, nover4);
        // a^e mod 4 = 3 if a = 3 mod 4 & e odd, 1 otherwise
        uint64_t resmod4 = 1 + 2*(((a & 3ULL) == 3ULL) && (e & 1ULL));
        // resmodnover4 + x*nover4 = resmod4 mod 4
        uint64_t x = (nover4*(resmod4 - resmodnover4)) & 3ULL;
        return x*nover4 + resmodnover4;
      } else {
        uint64_t nctz = __builtin_ctzll(n);
        uint64_t oddn = n >> nctz;
        uint64_t noveroddn = 1ULL << nctz;
        uint64_t mask = noveroddn-1ULL;
        uint64_t atoemododdn = modpowu64oddn(a, e, oddn);
        uint64_t atoemod2pow = modpowu64general(a & mask, (a&1ULL ? e % (noveroddn/2) : e), noveroddn);
#ifdef USEHACKERSDELIGHTFNS
        uint64_t invatoemododdn = modinv64x(oddn);
#else
        uint64_t invatoemododdn = modinv64(oddn);
#endif
        // a^e = atoemododdn + m*oddn mod n
        // m = (atoemod2pow - atoemododdn)*oddn^(-1) mod noveroddn
        uint64_t m = ((atoemod2pow - atoemododdn)*invatoemododdn) & mask;
        return atoemododdn + m*oddn;
      } 
    } 
  }
}
