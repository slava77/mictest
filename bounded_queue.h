#ifndef _bounded_queue_
#define _bounded_queue_

#include <queue>

template<class T, class Container = std::vector<T>,
         class Compare = std::less<typename Container::value_type>>
class bounded_queue : private std::priority_queue<T, Container, Compare> {
public:
  typedef std::priority_queue<T, Container, Compare> parent_type;
  typedef typename parent_type::container_type container_type;
  typedef typename parent_type::value_type value_type;
  typedef typename parent_type::reference reference;
  typedef typename parent_type::const_reference const_reference;
  typedef typename parent_type::size_type size_type;

  using parent_type::empty;
  using parent_type::size;
  using parent_type::top;
  using parent_type::emplace;
  using parent_type::pop;

private:
  using parent_type::push;

protected:
  using parent_type::c;
  using parent_type::comp;

public:
  template <class... Args> reference emplace_start(Args&&... args)
  {
    if (size() == c.capacity()) pop();
    c.emplace_back(std::forward<Args>(args)...); return c.back();
  }
  void emplace_finish() { std::push_heap(c.begin(), c.end(), comp); }
  void clear() { c.clear(); }
  void swap(container_type& x) { c.swap(x); }
  void reserve(size_type n) { c.reserve(n); }
  void maybe_push(const value_type& x) {
    if (size() < c.capacity()) {
      push(x);
    } else {
      if (x < top()) {
        pop();
        push(x);
      }
    }
  }
  Container& rep() { return c; }
};

#endif
