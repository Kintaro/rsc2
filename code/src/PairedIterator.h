#ifndef __PAIRED_ITERATOR_H_
#define __PAIRED_ITERATOR_H_

#include <algorithm>
#include <functional>
#include <iostream>

template <class Proxy> class ProxyIterator 
{
private:
  typedef ProxyIterator<Proxy> S;
  typedef typename Proxy::InnerIterator InnerIterator;

public:
  typedef std::random_access_iterator_tag iterator_category;
  typedef typename Proxy::value_type value_type;
  typedef std::ptrdiff_t difference_type;
  typedef Proxy & reference;
  typedef Proxy * pointer;

  ProxyIterator() {}

  // For cast from non const to const.
  template <class AlternateProxy> ProxyIterator(const ProxyIterator<AlternateProxy> &in) : p_(*in) {}
  explicit ProxyIterator(const Proxy &p) : p_(p) {}

  // p_'s swap does value swapping, but here we want iterator swapping
  friend inline void swap(ProxyIterator<Proxy> &first, ProxyIterator<Proxy> &second) 
  {
    swap(first.I(), second.I());
  }

  // p_'s operator= does value copying, but here we want iterator copying.
  S &operator=(const S &other) 
  {
    I() = other.I();
    return *this;
  }

  bool operator==(const S &other) const { return I() == other.I(); }
  bool operator!=(const S &other) const { return !(*this == other); }
  bool operator<(const S &other) const { return I() < other.I(); }
  bool operator>(const S &other) const { return other < *this; }
  bool operator<=(const S &other) const { return !(*this > other); }
  bool operator>=(const S &other) const { return !(*this < other); }

  S &operator++() { return *this += 1; }
  S operator++(int) { S ret(*this); ++*this; return ret; }
  S &operator+=(std::ptrdiff_t amount) { I() += amount; return *this; }
  S operator+(std::ptrdiff_t amount) const { S ret(*this); ret += amount; return ret; }

  S &operator--() { return *this -= 1; }
  S operator--(int) { S ret(*this); --*this; return ret; }
  S &operator-=(std::ptrdiff_t amount) { I() += (-amount); return *this; }
  S operator-(std::ptrdiff_t amount) const { S ret(*this); ret -= amount; return ret; }

  std::ptrdiff_t operator-(const S &other) const { return I() - other.I(); }

  Proxy &operator*() { return p_; }
  const Proxy &operator*() const { return p_; }
  Proxy *operator->() { return &p_; }
  const Proxy *operator->() const { return &p_; }
  Proxy operator[](std::ptrdiff_t amount) const { return *(*this + amount); }

  const InnerIterator &Inner() { return p_.Inner(); }

  private:
  InnerIterator &I() { return p_.Inner(); }
  const InnerIterator &I() const { return p_.Inner(); }

  Proxy p_;
};

template <class Proxy> ProxyIterator<Proxy> operator+(std::ptrdiff_t amount, const ProxyIterator<Proxy> &it) {
  return it + amount;
}

template <class KeyIter, class ValueIter> class JointProxy;

template <class KeyIter, class ValueIter> class JointIter 
{
public:
  JointIter() {}

  JointIter(const KeyIter &key_iter, const ValueIter &value_iter) : key_(key_iter), value_(value_iter) {}

  bool operator==(const JointIter<KeyIter, ValueIter> &other) const { return key_ == other.key_; }

  bool operator<(const JointIter<KeyIter, ValueIter> &other) const { return (key_ < other.key_); }

  std::ptrdiff_t operator-(const JointIter<KeyIter, ValueIter> &other) const { return key_ - other.key_; }

  JointIter<KeyIter, ValueIter> &operator+=(std::ptrdiff_t amount) 
  {
    key_ += amount;
    value_ += amount;
    return *this;
  }

  void swap(const JointIter &other) 
  {
    std::swap(key_, other.key_);
    std::swap(value_, other.value_);
  }

private:
  friend class JointProxy<KeyIter, ValueIter>;
  KeyIter key_;
  ValueIter value_;
};

template <class KeyIter, class ValueIter> class JointProxy 
{
private:
  typedef JointIter<KeyIter, ValueIter> InnerIterator;

public:
  typedef struct 
  {
    typename std::iterator_traits<KeyIter>::value_type key;
    typename std::iterator_traits<ValueIter>::value_type value;
    const typename std::iterator_traits<KeyIter>::value_type &GetKey() const { return key; }
  } value_type;

  JointProxy(const KeyIter &key_iter, const ValueIter &value_iter) : inner_(key_iter, value_iter) {}
  JointProxy(const JointProxy<KeyIter, ValueIter> &other) : inner_(other.inner_) {}

  operator value_type() const 
  {
    value_type ret;
    ret.key = *inner_.key_;
    ret.value = *inner_.value_;
    return ret;
  }

  JointProxy &operator=(const JointProxy &other) 
  {
    *inner_.key_ = *other.inner_.key_;
    *inner_.value_ = *other.inner_.value_;
    return *this;
  }

  JointProxy &operator=(const value_type &other) 
  {
    *inner_.key_ = other.key;
    *inner_.value_ = other.value;
    return *this;
  }

  typename std::iterator_traits<KeyIter>::reference GetKey() const 
  {
    return *(inner_.key_);
  }

  void swap(JointProxy<KeyIter, ValueIter> &other) 
  {
    std::swap(*inner_.key_, *other.inner_.key_);
    std::swap(*inner_.value_, *other.inner_.value_);
  }

private:
  friend class ProxyIterator<JointProxy<KeyIter, ValueIter> >;

  InnerIterator &Inner() { return inner_; }
  const InnerIterator &Inner() const { return inner_; }
  InnerIterator inner_;
};

template <class Proxy, class Less> class LessWrapper : public std::binary_function<const typename Proxy::value_type &, const typename Proxy::value_type &, bool> 
{
public:
  explicit LessWrapper(const Less &less) : less_(less) {}

  bool operator()(const Proxy &left, const Proxy &right) const 
  {
    return less_(left.GetKey(), right.GetKey());
  }
  bool operator()(const Proxy &left, const typename Proxy::value_type &right) const 
  {
    return less_(left.GetKey(), right.GetKey());
  }
  bool operator()(const typename Proxy::value_type &left, const Proxy &right) const 
  {
    return less_(left.GetKey(), right.GetKey());
  }
  bool operator()(const typename Proxy::value_type &left, const typename Proxy::value_type &right) const 
  {
    return less_(left.GetKey(), right.GetKey());
  }

private:
  const Less less_;
};

template <class KeyIter, class ValueIter> class PairedIterator : public ProxyIterator<JointProxy<KeyIter, ValueIter> > 
{
public:
  PairedIterator(const KeyIter &key, const ValueIter &value) :
  ProxyIterator<JointProxy<KeyIter, ValueIter> >(JointProxy<KeyIter, ValueIter>(key, value)) {}
};

template <class KeyIter, class ValueIter, class Less> void JointSort(const KeyIter &key_begin, const KeyIter &key_end, const ValueIter &value_begin, const Less &less) 
{
  ProxyIterator<JointProxy<KeyIter, ValueIter> > full_begin(JointProxy<KeyIter, ValueIter>(key_begin, value_begin));
  LessWrapper<JointProxy<KeyIter, ValueIter>, Less> less_wrap(less);
  std::sort(full_begin, full_begin + (key_end - key_begin), less_wrap);
}

template <class KeyIter, class ValueIter> void JointSort(const KeyIter &key_begin, const KeyIter &key_end, const ValueIter &value_begin) 
{
  JointSort(key_begin, key_end, value_begin, std::less<typename std::iterator_traits<KeyIter>::value_type>());
}

template <class KeyIter, class ValueIter, class ValueIterB, class Less> void JointSort3(const KeyIter &key_begin, const KeyIter &key_end, const ValueIter &value_begin, const ValueIterB &valueB_begin, const Less &less) 
{
  ProxyIterator<JointProxy<ValueIter, ValueIterB> > pair_begin(JointProxy<ValueIter, ValueIterB>(value_begin, valueB_begin));
  ProxyIterator<JointProxy<KeyIter, ProxyIterator<JointProxy<ValueIter, ValueIterB>>>> full_begin(JointProxy<KeyIter, ProxyIterator<JointProxy<ValueIter, ValueIterB> >>(key_begin, pair_begin));
  LessWrapper<JointProxy<KeyIter, ProxyIterator<JointProxy<ValueIter, ValueIterB> >>, Less> less_wrap(less);
  std::sort(full_begin, full_begin + (key_end - key_begin), less_wrap);
}

template <class KeyIter, class ValueIter, class ValueIterB> void JointSort3(const KeyIter &key_begin, const KeyIter &key_end, const ValueIter &value_begin, const ValueIterB &valueB_begin) 
{
  JointSort3(key_begin, key_end, value_begin, valueB_begin, std::less<typename std::iterator_traits<KeyIter>::value_type>());
}

template <class KeyIter, class ValueIter, class Less> void JointPartialSort(const KeyIter &key_begin, const KeyIter &key_middle, const KeyIter &key_end, const ValueIter &value_begin, const ValueIter &value_middle, const Less &less) 
{
  ProxyIterator<JointProxy<KeyIter, ValueIter> > full_begin(JointProxy<KeyIter, ValueIter>(key_begin, value_begin));
  ProxyIterator<JointProxy<KeyIter, ValueIter> > full_middle(JointProxy<KeyIter, ValueIter>(key_middle, value_middle));
  LessWrapper<JointProxy<KeyIter, ValueIter>, Less> less_wrap(less);
  std::partial_sort(full_begin, full_middle, full_begin + (key_end - key_begin), less_wrap);
}

template <class KeyIter, class ValueIter> void JointPartialSort(const KeyIter &key_begin, const KeyIter &key_end, const ValueIter &value_begin) 
{
  JointPartialSort(key_begin, key_end, value_begin, std::less<typename std::iterator_traits<KeyIter>::value_type>());
}

namespace std 
{
  template <class KeyIter, class ValueIter> void swap(JointIter<KeyIter, ValueIter> &left, JointIter<KeyIter, ValueIter> &right) 
  {
    left.swap(right);
  }

  template <class KeyIter, class ValueIter> void swap(JointProxy<KeyIter, ValueIter> &left, JointProxy<KeyIter, ValueIter> &right) 
  {
  left.swap(right);
  }
} 

#endif 
