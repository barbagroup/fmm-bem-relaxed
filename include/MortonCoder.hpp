#pragma once

#include "BoundingBox.hpp"

#include <assert.h>

/** @file MortonCoder.hpp
 * @brief Define the MortonCoder class for Z-order-curve values, aka Morton
 *   codes.
 */

/** @class MortonCoder
 * @brief Class representing Z-order-curve values, aka Morton codes.
 *
 * The Z-order curve is a space-filling curve: a one-dimensional curve that
 * fills a multi-dimensional space. Space-filling curves offer advantages for
 * representing points in 3D space. Points near each other in 3D space are
 * likely to have close Morton codes! So we can store points in a map with
 * Z-order value as key, then iterate over nearby Z-order values to fetch
 * points nearby in space.
 *
 * Unfortunately, it's impossible to reduce 3D space to 1D space perfectly:
 * there are some discontinuities in any space-filling curve mapping. But the
 * mapping is still an awesome tool, and some optimizations (see BIGMIN on the
 * Wikipedia page linked below) make them very effective.
 *
 * The MortonCoder class encapsulates a BoundingBox and can be used to translate
 * between spatial Points and Morton codes relative to that BoundingBox.
 *
 * A single Morton code corresponds to a rectangular volume within that
 * BoundingBox, called its <em>cell</em>. Each side of the BoundingBox is
 * divided into 2^L equal-sized cells, for a total of 8^L cells.
 *
 * This class computes maps box numbers to point and visa-versa
 * with respect to a bounding box and the number of equal-volume boxes (8^L).
 * These mappings are performed in O(1) time.
 */

// TODO: absorb into Octree...
template <typename Point>
class MortonCoder {
  // Using a 32-bit unsigned int for the code_type
  // means we can only resolve 10 3D levels
  static constexpr unsigned L = 10;
  static_assert(1 <= L && L <= 10, "L (LEVELS) must be between 1 and 10");

 public:

  typedef Point point_type;

  typedef unsigned code_type;

  /** The number of bits per dimension [octree subdivisions]. #cells = 8^L. */
  static constexpr unsigned levels = L;
  /** The number of cells per side of the bounding box (2^L). */
  static constexpr code_type cells_per_side = code_type(1) << L;
  /** One more than the largest code (8^L). */
  static constexpr code_type end_code = code_type(1) << (3*L);

  /** Construct a MortonCoder with a bounding box. */
  MortonCoder(const BoundingBox<point_type>& bb)
    : pmin_(bb.min()),
      cell_size_((bb.max() - bb.min()) / cells_per_side) {
    assert(!bb.empty());
  }

  /** Return the MortonCoder's bounding box. */
  BoundingBox<point_type> bounding_box() const {
    return BoundingBox<point_type>(pmin_, pmin_ + (cell_size_ * cells_per_side));
  }

  /** Return the bounding box of the cell with Morton code @a c.
   * @pre c < end_code */
  BoundingBox<point_type> cell(code_type c) const {
    assert(c < end_code);
    point_type p = deinterleave(c);
    p *= cell_size_;
    p += pmin_;
    return BoundingBox<point_type>(p, p + cell_size_);
  }

  /** Return the Morton code of Point @a p.
   * @pre bounding_box().contains(@a p)
   * @post cell(result).contains(@a p) */
  code_type code(const point_type& p) const {
    assert(bounding_box().contains(p));
    point_type s = (p - pmin_) / cell_size_;
    return interleave((unsigned) s[0], (unsigned) s[1], (unsigned) s[2]);
  }

  // 0x09249249 = 0b001001001001001001001001001001
  static constexpr code_type coordinate_mask = 0x09249249;
  static constexpr code_type x_mask = coordinate_mask << 0;
  static constexpr code_type y_mask = coordinate_mask << 1;
  static constexpr code_type z_mask = coordinate_mask << 2;

  /** True if min <= idx <= max and idx is inside the box defined
   * by the Morton codes @a min and @a max
   */
  bool is_in_box(code_type idx, code_type min, code_type max) const {
    return (min & x_mask) <= (idx & x_mask) && (idx & x_mask) <= (max & x_mask)
        && (min & y_mask) <= (idx & y_mask) && (idx & y_mask) <= (max & y_mask)
        && (min & z_mask) <= (idx & z_mask) && (idx & z_mask) <= (max & z_mask);
  }

  /** Advance idx to the next box contained in the bounding box defined
   * by the Morton codes min and max
   *
   * @return idx if is_in_box(idx,min,max)
   *         min if idx <= min
   *         idx if idx >= max
   *         new idx >= idx such that is_in_box(new idx,min,max)
   *                        and for no n in [idx,new idx), is_in_box(n,min,max)
   * @post result > max || is_in_box(result, min, max)
   * @post result >= idx
   */
  code_type advance_to_box(code_type idx, code_type min, code_type max) const {
    if (idx >= max) return idx;

    // If outside the box in some coord, record the difference
    code_type delta = 0;
    if      ((idx & x_mask) > (max & x_mask))  delta |= (idx ^ max) & x_mask;
    else if ((idx & x_mask) < (min & x_mask))  delta |= (idx ^ min) & x_mask;
    if      ((idx & y_mask) > (max & y_mask))  delta |= (idx ^ max) & y_mask;
    else if ((idx & y_mask) < (min & y_mask))  delta |= (idx ^ min) & y_mask;
    if      ((idx & z_mask) > (max & z_mask))  delta |= (idx ^ max) & z_mask;
    else if ((idx & z_mask) < (min & z_mask))  delta |= (idx ^ min) & z_mask;

    // Delta is only zero if idx is in the box
    if (delta == 0) return idx;

    // Smear into a low bit mask, i.e. 0000111111111111
    delta = smear_low_1(delta >> 1);
    if ((delta+1) & idx) {     // The idx bit above delta mask is one, need zero
      // Chi masks high bits we cannot carry into
      code_type chi = ~smear_low_3(idx ^ max);
      // The first 0 in idx and chi that is higher than delta
      delta = ~(idx | chi | delta);
      delta = (delta & -delta) - 1;
    }

    // Flip zero bit of idx and zero all lower bits
    idx = (idx | delta) + 1;

    // For each coordinate, if idx is low set to min
    if ((idx & x_mask) < (min & x_mask))  idx |= min & x_mask;
    if ((idx & y_mask) < (min & y_mask))  idx |= min & y_mask;
    if ((idx & z_mask) < (min & z_mask))  idx |= min & z_mask;

    return idx;
  }

 private:

  /** The minimum of the MortonCoder bounding box. */
  point_type pmin_;
  /** The extent of a single cell. */
  point_type cell_size_;

  /** Spreads the bits of a 10-bit number so that there are two 0s
   *  in between each bit.
   * @param x 10-bit integer
   * @return 28-bit integer of form 0b0000X00X00X00X00X00X00X00X00X00X,
   * where the X's are the original bits of @a x
   */
  inline unsigned spread_bits(unsigned x) const {
    x = (x | (x << 16)) & 0b00000011000000000000000011111111;
    x = (x | (x <<  8)) & 0b00000011000000001111000000001111;
    x = (x | (x <<  4)) & 0b00000011000011000011000011000011;
    x = (x | (x <<  2)) & 0b00001001001001001001001001001001;
    return x;
  }

  /** Interleave the bits of n into x, y, and z.
   * @pre x = [... x_2 x_1 x_0]
   * @pre y = [... y_2 y_1 y_0]
   * @pre z = [... z_2 z_1 z_0]
   * @post n = [... z_1 y_1 x_1 z_0 y_0 x_0]
   */
  inline code_type interleave(unsigned x, unsigned y, unsigned z) const {
    return spread_bits(x) | (spread_bits(y) << 1) | (spread_bits(z) << 2);
  }

  /** Does the inverse of spread_bits, extracting a 10-bit number from
   * a 28-bit number.
   * @param x 28-bit integer of form 0bYYYYXYYXYYXYYXYYXYYXYYXYYXYYXYYX
   * @return 10-bit integer of form 0b00...000XXXXXXXXXX,
   * where the X's are every third bit of @a x
   */
  inline unsigned compact_bits(unsigned x) const {
    x &= 0b00001001001001001001001001001001;
    x = (x | (x >>  2)) & 0b00000011000011000011000011000011;
    x = (x | (x >>  4)) & 0b00000011000000001111000000001111;
    x = (x | (x >>  8)) & 0b00000011000000000000000011111111;
    x = (x | (x >> 16)) & 0b00000000000000000000001111111111;
    return x;
  }

  /** Deinterleave the bits from n into a Point.
   * @pre n = [... n_2 n_1 n_0]
   * @post result.x = [... n_6 n_3 n_0]
   * @post result.y = [... n_7 n_4 n_1]
   * @post result.z = [... n_8 n_5 n_2]
   */
  inline point_type deinterleave(code_type c) const {
    typedef typename point_type::value_type value_type;
    return point_type(value_type(compact_bits(c)),
                      value_type(compact_bits(c >> 1)),
                      value_type(compact_bits(c >> 2)));
  }


  /** Smears the bits in c into the low bits by steps of one
   *
   * Example: 00011100100 -> 000111111111
   */
  inline code_type smear_low_1(code_type c) const {
    c |= c >>  1;
    c |= c >>  2;
    c |= c >>  4;
    c |= c >>  8;
    c |= c >> 16;
    return c;
  }

  /** Smears the bits in c into the low bits by steps of three
   *
   * Example: 0000010000000000 -> 0000010010010010
   */
  inline code_type smear_low_3(code_type c) const {
    c |= c >>  3;
    c |= c >>  6;
    c |= c >> 12;
    c |= c >> 24;
    return c;
  }
};
