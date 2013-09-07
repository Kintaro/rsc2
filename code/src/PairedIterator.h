#ifndef PAIRED_ITERATOR_H
#define PAIRED_ITERATOR_H

#include <iterator>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std;

/// @cond
template <typename ItA, typename ItB>
struct PairView {
    typedef typename ItA::value_type first_type;
    typedef typename ItB::value_type second_type;

    typedef std::pair<first_type, second_type> pair_type;

    PairView() {}
    PairView(const ItA &a, const ItB &b):
        first(*a), second(*b) {}

    PairView &operator=(const PairView &x)
        { first = x.first; second = x.second; return *this; }
    PairView &operator=(const pair_type &x)
        { first = x.first; second = x.second; return *this; }

    typename ItA::reference first;
    typename ItB::reference second;
    operator pair_type() const
        { return std::make_pair(first, second); }

    friend bool operator==(const PairView &a, const PairView &b)
        { return (a.first == b.first) && (a.second == b.second); }
    friend bool operator<(const PairView &a, const PairView &b)
        { return (a.first < b.first) || ((a.first == b.first) && (a.second < b.second)); }
    friend bool operator!=(const PairView &a, const PairView &b)
        { return !(a == b); }
    friend bool operator>(const PairView &a, const PairView &b)
        { return (b < a); }
    friend bool operator<=(const PairView &a, const PairView &b)
        { return !(b < a); }
    friend bool operator>=(const PairView &a, const PairView &b)
        { return !(a < b); }

    friend bool operator==(const PairView &a, const pair_type &b)
        { return (a.first == b.first) && (a.second == b.second); }
    friend bool operator<(const PairView &a, const pair_type &b)
        { return (a.first < b.first) || ((a.first == b.first) && (a.second < b.second)); }
    friend bool operator!=(const PairView &a, const pair_type &b)
        { return !(a == b); }
    friend bool operator>(const PairView &a, const pair_type &b)
        { return (b < a); }
    friend bool operator<=(const PairView &a, const pair_type &b)
        { return !(b < a); }
    friend bool operator>=(const PairView &a, const pair_type &b)
        { return !(a < b); }

    friend bool operator==(const pair_type &a, const pair_type &b)
        { return (a.first == b.first) && (a.second == b.second); }
    friend bool operator<(const pair_type &a, const pair_type &b)
        { return (a.first < b.first) || ((a.first == b.first) && (a.second < b.second)); }
    friend bool operator!=(const pair_type &a, const pair_type &b)
        { return !(a == b); }
    friend bool operator>(const pair_type &a, const pair_type &b)
        { return (b < a); }
    friend bool operator<=(const pair_type &a, const pair_type &b)
        { return !(b < a); }
    friend bool operator>=(const pair_type &a, const pair_type &b)
        { return !(a < b); }
};

/*!
 * @internal
 */
template <typename ItA, typename ItB>
struct PairedIterator {
    // --- standard iterator traits
    typedef typename PairView<ItA, ItB>::pair_type value_type;
    typedef PairView<ItA, ItB> reference;
    typedef PairedIterator<ItA,ItB> pointer;

    typedef typename std::iterator_traits<ItA>::difference_type difference_type;
    typedef std::random_access_iterator_tag iterator_category;

    // --- methods not required by the Random Access Iterator concept
    PairedIterator(const ItA &a, const ItB &b):
        a(a), b(b) {}

    // --- iterator requirements

    // default construction
    PairedIterator() {}

    // copy construction and assignment
    PairedIterator(const PairedIterator &x):
        a(x.a), b(x.b) {}
    PairedIterator &operator=(const PairedIterator &x)
        { a = x.a; b = x.b; return *this; }

    // pre- and post-increment
    PairedIterator &operator++()
        { ++a; ++b; return *this; }
    PairedIterator operator++(int)
        { PairedIterator tmp(*this); ++(*this); return tmp; }

    // pre- and post-decrement
    PairedIterator &operator--()
        { --a; --b; return *this; }
    PairedIterator operator--(int)
        { PairedIterator tmp(*this); --(*this); return tmp; }

    // arithmetic
    PairedIterator &operator+=(const difference_type &n)
        { a += n; b += n; return *this; }
    friend PairedIterator operator+(const PairedIterator &x, const difference_type &n)
        { return PairedIterator(x.a+n, x.b+n); }
    friend PairedIterator operator+(const difference_type &n, const PairedIterator &x)
        { return PairedIterator(x.a+n, x.b+n); }
    PairedIterator &operator-=(const difference_type &n)
        { a -= n; b -= n; return *this; }
    friend PairedIterator operator-(const PairedIterator &x, const difference_type &n)
        { return PairedIterator(x.a-n, x.b-n); }
    friend difference_type operator-(const PairedIterator &x, const PairedIterator &y)
        { return (x.a - y.a); }

    // (in-)equality and ordering
    friend bool operator==(const PairedIterator &x, const PairedIterator &y)
        { return (x.a == y.a) && (x.b == y.b); }
    friend bool operator<(const PairedIterator &x, const PairedIterator &y)
        { return (x.a < y.a); }

    // derived (in-)equality and ordering operators
    friend bool operator!=(const PairedIterator &x, const PairedIterator &y)
        { return !(x == y); }
    friend bool operator>(const PairedIterator &x, const PairedIterator &y)
        { return (y < x); }
    friend bool operator<=(const PairedIterator &x, const PairedIterator &y)
        { return !(y < x); }
    friend bool operator>=(const PairedIterator &x, const PairedIterator &y)
        { return !(x < y); }

    // dereferencing and random access
    reference operator*() const
        { return reference(a,b); }
    reference operator[](const difference_type &n) const
        { return reference(a+n, b+n); }
private:
    ItA a;
    ItB b;
};

template <typename ItA, typename ItB>
PairedIterator<ItA, ItB> make_paired_iterator(const ItA &a, const ItB &b)
{ return PairedIterator<ItA, ItB>(a, b); }
/// @endcond

#endif /* PAIRED_ITERATOR_H */