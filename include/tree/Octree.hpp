#pragma once

#include <iostream>
#include <iomanip>
#include <assert.h>

#include "BoundingBox.hpp"
#include "Util.hpp"

// Automatically derive !=, <=, >, and >= from a class's == and <
using namespace std::rel_ops;


//! Class for tree structure
template <typename Source, typename Point>
class Octree
{
 public:
  // Type declarations
  typedef Source source_type;
  typedef Point point_type;
  static_assert(point_type::dimension == 3, "Only 3D at the moment");

  //! The type of this tree
  typedef Octree<source_type, point_type> tree_type;

 private:
  template <typename SourceIter>
  BoundingBox<point_type> get_boundingbox(SourceIter begin, SourceIter end) {
    BoundingBox<point_type> result;
    for ( ; begin != end; ++begin)
      result |= static_cast<point_type>(*begin);
    // Make sure the bounding box is square and slightly scaled
    // TODO: improve
    auto dim = result.dimensions();
    auto maxdim = std::max(dim[0], std::max(dim[1], dim[2]));
    result |= result.min() + point_type(maxdim) * (1 + 1e-6);
    //std::cout << "Bounding Box: " << result << "\n";
    return result;
  }

  // The Coder this tree is based on
  class MortonCoder {
   public:
    // Using a 32-bit unsigned int for the code_type
    // means we can only resolve 10 3D levels
    typedef unsigned code_type;

    /** The number of bits per dimension [octree subdivisions]. #cells = 8^L. */
    static constexpr unsigned levels = 10;
    /** The number of cells per side of the bounding box (2^L). */
    static constexpr code_type cells_per_side = code_type(1) << levels;
    /** One more than the largest code (8^L). */
    static constexpr code_type end_code = code_type(1) << (3*levels);

    /** Construct a MortonCoder with a bounding box. */
    MortonCoder(const BoundingBox<point_type>& bb)
        : pmin_(bb.min()),
          cell_size_((bb.max() - bb.min()) / cells_per_side) {
      assert(!bb.empty());
    }

    MortonCoder() {};

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
      point_type s = (p - pmin_) / cell_size_;
      assert((unsigned) s[0] < cells_per_side &&
             (unsigned) s[1] < cells_per_side &&
             (unsigned) s[2] < cells_per_side);
      return interleave((unsigned) s[0], (unsigned) s[1], (unsigned) s[2]);
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
  };

  // Code type
  typedef typename MortonCoder::code_type code_type;

  // Morton coder to use for the points
  MortonCoder coder_;

  // Morton coded objects this Tree holds.
  std::vector<source_type> source_;
  std::vector<point_type> point_;
  std::vector<code_type> mc_;
  std::vector<unsigned> permute_;
  std::vector<unsigned> level_offset_;

  struct box_data {
    static constexpr unsigned leaf_bit = (1<<31);
    static constexpr unsigned max_marker_bit = (1<<30);

    //! key_ = leaf_bit 0* marker_bit morton_code
    unsigned key_;   // TODO: use level + leaf_bit instead of key
    unsigned parent_;
    // These can be either point offsets or box offsets depending on is_leaf
    unsigned child_begin_;
    unsigned child_end_;
    // TODO: body_begin_ and body_end_?

    box_data(unsigned key, unsigned parent=0,
             unsigned child_begin=0, unsigned child_end=0)
        : key_(key), parent_(parent),
          child_begin_(child_begin), child_end_(child_end) {
    }

    unsigned num_children() const {
      return child_end_ - child_begin_;
    }

    /** Gets the level from the key
     * TODO: optimize
     */
    unsigned level() const {
      static constexpr unsigned lookup[] = {0, 3, 0, 3, 4, 7, 0, 9,
                                            3, 4, 5, 6, 7, 8, 1, 10,
                                            2, 4, 6, 9, 5, 5, 8, 2,
                                            6, 9, 7, 2, 8, 1, 1, 10};
      unsigned v = key_ & ~leaf_bit;
      v |= v >> 1;
      v |= v >> 2;
      v |= v >> 4;
      v |= v >> 8;
      v |= v >> 16;
      return lookup[(v * 0X07C4ACDD) >> 27];
    }

    /** Returns the minimum possible Morton code in this box
     * TODO: this is a stupid way of doing this
     */
    code_type get_mc_lower_bound() const {
      code_type mc_mask = key_;
      while (!(mc_mask & max_marker_bit))
        mc_mask = mc_mask << 3;
      return mc_mask & ~max_marker_bit;
    }
    /** Returns the maximum possible Morton code in this box
     * TODO: this is a stupid way of doing this
     */
    code_type get_mc_upper_bound() const {
      code_type mc_mask = key_;
      while (!(mc_mask & max_marker_bit))
        mc_mask = (mc_mask << 3) | 7;
      return mc_mask & ~max_marker_bit;
    }

    void set_leaf(bool b) {
      if (b)
        key_ |= leaf_bit;
      else
        key_ &= ~leaf_bit;
    }

    bool is_leaf() const {
      return key_ & leaf_bit;
    }
  };

  std::vector<box_data> box_data_;

 public:
  // Predeclarations
  struct Body;
  typedef Body body_type;
  struct Box;
  typedef Box box_type;
  struct body_iterator;
  struct box_iterator;

  struct Body {
    /** Construct an invalid Body */
    Body() : idx_(0), tree_(NULL) {
    }

    const source_type& source() const {
      return tree_->source_[idx_];
    }
    source_type& source() {
      return tree_->source_[idx_];
    }
    const point_type& point() const {
      return tree_->point_[idx_];
    }
    point_type& point() {
      return tree_->point_[idx_];
    }
    //! The original order this body was seen
    unsigned number() const {
      return tree_->permute_[idx_];
    }
    //! A body's index in the tree
    unsigned index() const {
      return idx_;
    }
    code_type morton_index() const {
      return tree_->mc_[idx_];
    }

   private:
    unsigned idx_;
    tree_type* tree_;
    Body(unsigned idx, tree_type* tree)
        : idx_(idx), tree_(tree) {
      assert(idx_ < tree_->size());
    }
    friend class Octree;
  };

  // A tree-aligned box
  struct Box {
    /** Construct an invalid Box */
    Box() : idx_(0), tree_(NULL) {
    }

    unsigned index() const {
      return idx_;
    }
    code_type morton_index() const {
      return data().key_;
    }
    unsigned level() const {
      return data().level();
    }
    double side_length() const {
      return tree_->coder_.bounding_box().dimensions()[0] / (1 << level());
    }
    double radius() const {
      return side_length() / 2.0;
    }
    unsigned num_children() const {
      return data().num_children();
    }
    bool is_leaf() const {
      return data().is_leaf();
    }
    // TODO: optimize
    point_type center() const {
      BoundingBox<point_type> bb = tree_->coder_.cell(data().get_mc_lower_bound());
      point_type p = bb.min();
      p += bb.dimensions() * (1 << (10-data().level()-1));
      return p;
    }

    /** The parent box of this box */
    Box parent() const {
      return Box(data().parent_, tree_);
    }

    /** The begin iterator to the Points contained in this box */
    body_iterator body_begin() const {
      if (is_leaf())
        return body_iterator(data().child_begin_, tree_);
      else
	return child_begin()->body_begin();
    }
    /** The end iterator to the Points contained in this box */
    body_iterator body_end() const {
      if (is_leaf())
        return body_iterator(data().child_end_, tree_);
      else
	return (--child_end())->body_end();
    }

    /** The begin iterator to the child Boxes contained in this box */
    box_iterator child_begin() const {
      assert(!is_leaf());
      return box_iterator(data().child_begin_, tree_);
    }
    /** The end iterator to the child Boxes contained in this box */
    box_iterator child_end() const {
      assert(!is_leaf());
      return box_iterator(data().child_end_, tree_);
    }

    /** Write a Box to an output stream */
    inline friend std::ostream& operator<<(std::ostream& s,
					   const box_type& b) {
      return s << "Box " << b.index()
	       << " (Level " << b.level() << ", Parent " << b.parent().index()
	       << ", Bodies " << b.body_begin()->index()
	       << "-" << (--b.body_end())->index()
	       << "): " << b.center();
    }
   private:
    unsigned idx_;
    tree_type* tree_;
    Box(unsigned idx, tree_type* tree)
        : idx_(idx), tree_(tree) {
    }
    inline box_data& data() const {
      return tree_->box_data_[idx_];
    }
    friend class Octree;
  };

  /** @struct Tree::box_iterator
   * @brief Iterator class for Boxes in the tree
   * TODO: Use a Mutator to condense/clarify code
   */
  struct box_iterator {
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Box value_type;
    /** Type of pointers to elements. */
    typedef Box* pointer;
    /** Type of references to elements. */
    typedef Box& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Iterator difference */
    typedef std::ptrdiff_t difference_type;
    /** Construct an invalid box_iterator */
    box_iterator() : idx_(0), tree_(NULL) {
    }

    box_iterator& operator--() {
      --idx_;
      return *this;
    }
    box_iterator& operator++() {
      ++idx_;
      return *this;
    }
    box_iterator& operator+(int n) {
      idx_ += n;
      return *this;
    }
    box_iterator& operator-(int n) {
      idx_ -= n;
      return *this;
    }
    Box operator*() const {
      return Box(idx_, tree_);
    }
    Box* operator->() const {
      placeholder_ = operator*();
      return &placeholder_;
    }
    bool operator==(const box_iterator& it) const {
      return tree_ == it.tree_ && idx_ == it.idx_;
    }

   private:
    unsigned idx_;
    tree_type* tree_;
    mutable Box placeholder_;
    box_iterator(unsigned idx, tree_type* tree)
        : idx_(idx), tree_(tree) {
    }
    box_iterator(Box b)
        : idx_(b.idx_), tree_(b.tree_) {
    }
    friend class Octree;
  };

  /** @struct Tree::body_iterator
   * @brief Iterator class for Bodies in the tree
   * TODO: Use a Mutator to condense/clarify code
   */
  struct body_iterator {
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Body value_type;
    /** Type of pointers to elements. */
    typedef Body* pointer;
    /** Type of references to elements. */
    typedef Body& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Iterator difference */
    typedef std::ptrdiff_t difference_type;
    /** Construct an invalid iterator */
    body_iterator() : idx_(0), tree_(0) {
    }

    body_iterator& operator++() {
      ++idx_;
      return *this;
    }
    body_iterator& operator--() {
      --idx_;
      return *this;
    }
    body_iterator& operator+(int n) {
      idx_ += n;
      return *this;
    }
    body_iterator& operator-(int n) {
      idx_ -= n;
      return *this;
    }
    Body operator*() const {
      return Body(idx_, tree_);
    }
    Body* operator->() const {
      placeholder_ = operator*();
      return &placeholder_;
    }
    bool operator==(const body_iterator& it) const {
      return tree_ == it.tree_ && idx_ == it.idx_;
    }

   private:
    unsigned idx_;
    tree_type* tree_;
    mutable Body placeholder_;
    body_iterator(unsigned idx, tree_type* tree)
        :idx_(idx), tree_(tree) {
    }
    friend class Octree;
  };

  //! Construct an octree encompassing a bounding box
  Octree(const BoundingBox<Point>& bb)
    : coder_(bb) {
  }

  template <typename PointIter, typename Options>
  Octree(PointIter first, PointIter last,
	 Options& opts)
    : coder_(get_boundingbox(first, last)) {
    construct_tree(first, last, opts.max_per_box());
  }

  Octree() {};

  /** Return the Bounding Box that this Octree encompasses
   */
  BoundingBox<point_type> bounding_box() const {
    return coder_.bounding_box();
  }

  /** The number of points contained in this tree
   */
  inline unsigned size() const {
    return point_.size();
  }

  /** The number of points contained in this tree
   */
  inline unsigned bodies() const {
    return size();
  }

  /** The number of boxes contained in this tree
   */
  inline unsigned boxes() const {
    return box_data_.size();
  }

  /** The maximum level of any box in this tree
   */
  inline unsigned levels() const {
    return level_offset_.size() - 1;
  }

#if 0
  template <typename PointIter>
  void construct_tree(PointIter p_begin, PointIter p_end, unsigned NCRIT = 126) {
    // Create a code-idx pair vector
    typedef std::pair<code_type, unsigned> code_pair;
    std::vector<code_pair> codes;
    unsigned idx = 0;
    for (PointIter pi = p_begin; pi != p_end; ++pi, ++idx) {
      assert(coder_.bounding_box().contains(*pi));
      codes.push_back(std::make_pair(coder_.code(*pi), idx));
    }

    // TODO: Use radix sort for efficiency or incrementally sort...
    std::sort(codes.begin(), codes.end());

    std::vector<point_type> points_tmp(p_begin, p_end);
    // Extract the code, permutation vector, and sorted point
    for (auto it = codes.begin(); it != codes.end(); ++it) {
      mc_.push_back(it->first);
      permute_.push_back(it->second);
      point_.push_back(points_tmp[permute_.back()]);
    }

    // Push the root box which contains all points
    box_data_.push_back(box_data(1, 0, 0, mc_.size()));
    level_offset_.push_back(0);

    // For every box that is created
    // TODO: Can do this in one scan through the morton codes...
    for (unsigned k = 0; k != box_data_.size(); ++k) {
      if (box_data_[k].num_children() <= NCRIT) {
        box_data_[k].set_leaf(true);
        continue;
      }

      // Get the key and interval
      auto mc_begin = mc_.begin() + box_data_[k].child_begin_;
      auto mc_end   = mc_.begin() + box_data_[k].child_end_;
      unsigned shift    = 3*(MortonCoder::levels - box_data_[k].level() - 1);
      unsigned box_code = (*mc_begin) & (code_type(-1) << (shift+1));

      // Find the child box offsets
      // Construct off such that off[i],off[i+1] are the begin,end of child i
      std::vector<typename std::vector<code_type>::iterator> off(9);
      off[0] = mc_begin;
      off[8] = mc_end;
      off[4] = std::upper_bound(off[0], off[8], box_code | (4 << shift));
      off[2] = std::upper_bound(off[0], off[4], box_code | (2 << shift));
      off[1] = std::upper_bound(off[0], off[2], box_code | (1 << shift));
      off[3] = std::upper_bound(off[2], off[4], box_code | (3 << shift));
      off[6] = std::upper_bound(off[4], off[8], box_code | (6 << shift));
      off[5] = std::upper_bound(off[4], off[6], box_code | (5 << shift));
      off[7] = std::upper_bound(off[6], off[8], box_code | (7 << shift));

      // Split this box - point offsets become box offsets
      box_data_[k].child_begin_ = box_data_.size();
      box_data_[k].child_end_   = box_data_.size();

      // For each bucket
      for (int c = 0; c < 8; ++c) {
        unsigned begin_c = off[c]   - mc_.begin();
        unsigned end_c   = off[c+1] - mc_.begin();

        // If this child contains points, add this child box
        if (end_c - begin_c > 0) {
          // Construct the new box
          code_type key_c = (box_data_[k].key_ << 3) | c;
          box_data box_c(key_c, k, begin_c, end_c);

          // TODO: Optimize on key
          // If this is starting a new level, record it
          if (box_c.level() > levels())
            level_offset_.push_back(box_data_.size());

          // Increment parent child offset
          ++box_data_[k].child_end_;
          // Add the child
          box_data_.push_back(box_c);
        }
      }
    }

    level_offset_.push_back(box_data_.size());
  }
#endif

#if 1
  template <typename SourceIter>
  void construct_tree(SourceIter p_begin, SourceIter p_end, unsigned NCRIT = 126) {
    // create a new vector with all points instead of sources
    std::vector<point_type> source_points(p_end-p_begin);
    std::transform(p_begin, p_end, source_points.begin(), [](source_type si) { return static_cast<point_type>(si); });

    // Create a code-idx pair vector
    typedef std::pair<code_type, unsigned> code_pair;
    std::vector<code_pair> codes;
    unsigned idx = 0;
    for (auto pi = source_points.begin(); pi != source_points.end(); ++pi, ++idx) {
      assert(coder_.bounding_box().contains(*pi));
      codes.push_back(std::make_pair(coder_.code(*pi), idx));
    }

    // Push the root box which contains all points
    box_data_.push_back(box_data(1, 0, 0, codes.size()));
    level_offset_.push_back(0);

    // For every box that is created
    for (unsigned k = 0; k != box_data_.size(); ++k) {

      // If this box is has few enough points, mark as leaf and continue
      if (box_data_[k].num_children() <= NCRIT) {
        box_data_[k].set_leaf(true);
        continue;
      }

      // Get the box data
      auto code_begin = codes.begin() + box_data_[k].child_begin_;
      auto code_end   = codes.begin() + box_data_[k].child_end_;
      unsigned shift  = 3*(MortonCoder::levels - box_data_[k].level() - 1);

      // Sort the points in this box into the eight "bucket" children
      auto off = bucket_sort(code_begin, code_end, 8,
                             [shift] (code_pair& v)
                             { return (v.first >> shift) & 7; });

      // Split this box - point offsets become box offsets
      box_data_[k].child_begin_ = box_data_.size();
      box_data_[k].child_end_   = box_data_.size();

      // For each bucket
      for (int c = 0; c < 8; ++c) {
        unsigned begin_c = off[c]   - codes.begin();
        unsigned end_c   = off[c+1] - codes.begin();

        // If this child contains points, add this child box
        if (end_c - begin_c > 0) {
          // Construct the new box
          code_type key_c = (box_data_[k].key_ << 3) | c;
          box_data box_c(key_c, k, begin_c, end_c);

          // TODO: Optimize on key
          // If this is starting a new level, record it
          if (box_c.level() > levels())
            level_offset_.push_back(box_data_.size());

          // Increment parent child offset
          ++box_data_[k].child_end_;
          // Add the child
          box_data_.push_back(box_c);
        }
      }
    }

    level_offset_.push_back(box_data_.size());

    // Copy the points to a vector
    std::vector<point_type> points_tmp(source_points.begin(), source_points.end());
    std::vector<source_type> sources_tmp(p_begin, p_end);
    // Extract the code, permutation vector, and sorted point
    for (auto it = codes.begin(); it != codes.end(); ++it) {
      mc_.push_back(it->first);
      permute_.push_back(it->second);
      point_.push_back(points_tmp[permute_.back()]);
      source_.push_back(sources_tmp[permute_.back()]);
    }
  }
#endif

  /** Return the root box of this tree */
  Box root() const {
    return Box(0, const_cast<tree_type*>(this));
  }
  /** Return an iterator to the first body in this tree */
  body_iterator body_begin() const {
    return body_iterator(0, const_cast<tree_type*>(this));
  }
  /** Return an iterator one past the last body in this tree */
  body_iterator body_end() const {
    return body_iterator(point_.size(), const_cast<tree_type*>(this));
  }
  /** Return an iterator to the first box in this tree */
  box_iterator box_begin() const {
    return box_iterator(0, const_cast<tree_type*>(this));
  }
  /** Return an iterator one past the last box in this tree */
  box_iterator box_end() const {
    return box_iterator(box_data_.size(), const_cast<tree_type*>(this));
  }
  /** Return an iterator to the first box at level L in this tree
   * @pre L < levels()
   */
  box_iterator box_begin(unsigned L) const {
    assert(L < levels());
    return box_iterator(level_offset_[L], const_cast<tree_type*>(this));
  }
  /** Return an iterator one past the last box at level L in this tree
   * @pre L < levels()
   */
  box_iterator box_end(unsigned L) const {
    assert(L < levels());
    return box_iterator(level_offset_[L+1], const_cast<tree_type*>(this));
  }

  /** Permute a vector to the same order of the input points.
   *
   * @param[in] v The vector associated with the original input points
   * @returns A vector whose elements have been permuted into the same
   * order as the points in this tree.
   */
  template <typename T>
  std::vector<T> permute(const std::vector<T>& v) {
    std::vector<T> temp(v.size());
    for (unsigned i=0; i < v.size(); ++i)
      temp[i] = v[permute_[i]];
    return temp;
  }

  /** Inverse permute a vector from the same order of the input points.
   *
   * @param[in] v The vector associated with the current points in the tree
   * @returns A vector whose elements have been permuted into the same
   * order as the original input points of this tree.
   */
  template <typename T>
  std::vector<T> ipermute(const std::vector<T>& v) {
    std::vector<T> temp(v.size());
    for (unsigned i = 0; i < v.size(); ++i)
      temp[permute_[i]] = v[i];
    return temp;
  }

  /** Write an Octree to an output stream */
  inline friend std::ostream& operator<<(std::ostream& s,
					 const tree_type& t) {
    struct {
      inline std::ostream& print(std::ostream& ss,
				 const box_type& b) {
	ss << std::string(2*b.level(), ' ') << b;
	if (!b.is_leaf()) {
	  for (auto ci = b.child_begin(); ci != b.child_end(); ++ci) {
	    ss << "\n";
	    print(ss,*ci);
	  }
	}
	return ss;
      }
    } level_traverse;

    return level_traverse.print(s, t.root());
  }
};





