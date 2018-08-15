#ifndef RAN_H
#define RAN_H

#include "nr3.h"

struct Ran
{
  Ullong u, v, w;
  Ran(Ullong j) : v(4101842887655102017LL), w(1)
  {
    u = j ^ v;
    int64();
    v = u;
    int64();
    w = v;
    int64();
  }
  inline Ullong int64()
  {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17;
    v ^= v << 31;
    v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    Ullong x = u ^ (u << 21);
    x ^= x >> 35;
    x ^= x << 4;
    return (x + v) ^ w;
  }
  inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
  inline Uint int32() { return (Uint)int64(); }
};
struct Ranq1
{
  Ullong v;
  Ranq1(Ullong j) : v(4101842887655102017LL)
  {
    v ^= j;
    v = int64();
  }
  inline Ullong int64()
  {
    v ^= v >> 21;
    v ^= v << 35;
    v ^= v >> 4;
    return v * 2685821657736338717LL;
  }
  inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
  inline Uint int32() { return (Uint)int64(); }
};

struct Ranq2
{
  Ullong v, w;
  Ranq2(Ullong j) : v(4101842887655102017LL), w(1)
  {
    v ^= j;
    w = int64();
    v = int64();
  }
  inline Ullong int64()
  {
    v ^= v >> 17;
    v ^= v << 31;
    v ^= v >> 8;
    w = 4294957665U * (w & 0xffffffff) + (w >> 32);
    return v ^ w;
  }
  inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
  inline Uint int32() { return (Uint)int64(); }
};

struct Normaldev : Ran
{
  Doub mu, sig;
  Normaldev(Doub mmu, Doub ssig, Ullong i)
      : Ran(i), mu(mmu), sig(ssig) {}
  Doub dev()
  {
    Doub u, v, x, y, q;
    do
    {
      u = doub();
      v = 1.7156 * (doub() - 0.5);
      x = u - 0.449871;
      y = abs(v) + 0.386595;
      q = SQR(x) + y * (0.19600 * y - 0.25472 * x);
    } while (q > 0.27597 && (q > 0.27846 || SQR(v) > -4. * log(u) * SQR(u)));
    return mu + sig * v / u;
  }
};

#endif //RAN_H