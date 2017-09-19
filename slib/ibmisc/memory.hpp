/*
 * IBMisc: Misc. Routines for IceBin (and other code)
 * Copyright (c) 2013-2016 by Elizabeth Fischer
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// See: http://stackoverflow.com/questions/11002641/dynamic-casting-for-unique-ptr
#include <functional>
#include <memory>
#include <typeinfo>
#include <boost/variant.hpp>

namespace ibmisc {

/** Create a unique_ptr from a function that returns a value.  Eg:
    auto xptr(new_unique_ptr(make_x(...)));
*/

#if 0
template<class ValT>
std::unique_ptr<ValT> new_unique_ptr(ValT &&val)
    { return std::unique_ptr<ValT>(new ValT(std::move(val))); }

template<class ValT, typename... Args>
std::unique_ptr<ValT> new_unique_ptr(Args... args)
    { return std::unique_ptr<ValT>(new ValT(args...)); }
#endif


/** Convenience function to store an rvalue reference in a shared_ptr */
template<class ValT>
void reset_ptr(std::unique_ptr<ValT> &ptr, ValT &&val)
{ ptr.reset(new ValT(std::move(val))); }


template<class ValT>
std::unique_ptr<ValT> new_unique_ptr(ValT &&val)
    { return std::unique_ptr<ValT>(new ValT(std::move(val))); }


template<class ValT>
std::shared_ptr<ValT> new_shared_ptr()
    { return std::shared_ptr<ValT>(new ValT); }

// --------------------------------------------------------
// http://ficksworkshop.com/blog/14-coding/86-how-to-static-cast-std-unique-ptr 
template<typename D, typename B>
/** Static-cast a unique_ptr while moving it to another unique_ptr */
void static_move(std::unique_ptr<D> &dest, std::unique_ptr<B>& base)
{
    dest.reset(static_cast<D*>(base.release()));
}

template<typename D, typename B>
std::unique_ptr<D> static_cast_unique_ptr(std::unique_ptr<B>& base)
{
    return std::unique_ptr<D>(static_cast<D*>(base.release()));
}
  
template<typename D, typename B>
std::unique_ptr<D> static_cast_unique_ptr(std::unique_ptr<B>&& base)
{
    return std::unique_ptr<D>(static_cast<D*>(base.release()));
}
// --------------------------------------------------------
/** Dynamic-cast a unique_ptr while moving it to another unique_ptr */
template<typename D, typename B>
void dynamic_move(std::unique_ptr<D> &dest, std::unique_ptr<B>& base)
{
    dest.reset(dynamic_cast<D*>(base.release()));
}

template<typename D, typename B>
std::unique_ptr<D> dynamic_cast_unique_ptr(std::unique_ptr<B>& base)
{
    return std::unique_ptr<D>(dynamic_cast<D*>(base.release()));
}
  
template<typename D, typename B>
std::unique_ptr<D> dynamic_cast_unique_ptr(std::unique_ptr<B>&& base)
{
    return std::unique_ptr<D>(dynamic_cast<D*>(base.release()));
}
// --------------------------------------------------------



template <class T_DEST, class T_SRC>
inline std::shared_ptr<T_DEST> dynamic_shared_cast(std::unique_ptr<T_SRC> &&src)
{
    return std::shared_ptr<T_DEST>(dynamic_cast_unique_ptr<T_DEST, T_SRC>(std::move(src)));
}

// ---------------------------------------------------------------

/** Stores a value that is lazily evaluted.  The evaluation function
may return either a object (and "owned" reference), or a C-pointer (a
"borrowed" reference). */
template<class TypeT>
class LazyPtr {
    TypeT *_ptr;                        // Borrowed reference
    std::unique_ptr<TypeT> _uptr;   // Owned reference

    boost::variant<
        std::function<TypeT *()>,       // Borrowed
        std::function<std::unique_ptr<TypeT> ()>    // Owned
    > _compute;


public:

    /** Constructs with an already-evaluted borrowed reference. */
    LazyPtr(TypeT *ptr) : _ptr(ptr) {}

    /** Constructs with an already-evaluated owned refernece. */
    LazyPtr(std::unique_ptr<TypeT> &&uptr) : _uptr(std::move(uptr)) {
        _ptr = _uptr.get();
    }

    /** Constructs with a function to produce an owned reference. */
    LazyPtr(std::function<std::unique_ptr<TypeT> ()> const &compute_owned)
        : _ptr(0), _compute(compute_owned) {}

    /** Constructs with a function to produce a borrowed referene. */
    LazyPtr(std::function<TypeT *()> const &compute_borrowed)
        : _ptr(0), _compute(compute_borrowed) {}




    /** Constructs with a function to produce an owned reference. */
    LazyPtr(std::function<std::unique_ptr<TypeT> ()> &&compute_owned)
        : _ptr(0), _compute(std::move(compute_owned)) {}

    /** Constructs with a function to produce a borrowed referene. */
    LazyPtr(std::function<TypeT *()> &&compute_borrowed)
        : _ptr(0), _compute(std::move(compute_borrowed)) {}

    // -----------------------------------------------------
    // Used by operator*()
    class eval_visitor : public boost::static_visitor<> {
        LazyPtr<TypeT> const *lptr;
    public:
        eval_visitor(LazyPtr<TypeT> const *_lptr) : lptr(_lptr) {}

        void operator()(std::function<TypeT *()> const &borrowed_fn) const
        {
            const_cast<TypeT *&>(lptr->_ptr) = borrowed_fn(); }
        void operator()(std::function<std::unique_ptr<TypeT> ()> const &owned_fn) const {
            const_cast<std::unique_ptr<TypeT>&>(lptr->_uptr) = owned_fn();
            const_cast<TypeT *&>(lptr->_ptr) = &*lptr->_uptr;
        }
    };
    friend class eval_visitor;


    /** Dereference our "pointer."  Evaluates the function, if it has
    not already been evaluated. */
    TypeT &operator*() const {
        if (!_ptr) {
            boost::apply_visitor(eval_visitor(this), _compute);
        }
        return *_ptr;
    }
    TypeT *operator->() const
        { return &operator*(); }
};
// --------------------------------------------------------------

/** An allocatoer that retains ownership of everything it allocates;
and de-allocates all those things when it is destroyed. */
class TmpAlloc {
    std::vector<std::function<void()>> deleters;

    template<class T>
    static void del(T *t)
        { delete t; }

public:
    template<class T, typename... Args>
    T *newptr(Args... args)
    {
        T *ret = new T(args...);
        deleters.push_back(std::bind(&TmpAlloc::del<T>, ret));
        return ret;
    }

	/** Allocates, and stores a pointer to it in ptr.
	This template is eaiser to use than newptr() above, requires no type tags */
    template<class T, typename... Args>
    void reset_ptr(T *&ptr, Args... args)
    {
        ptr = new T(args...);
        deleters.push_back(std::bind(&TmpAlloc::del<T>, ptr));
    }


    /** Allocate and run an arbitrary constructor */
    template<class T, typename... Args>
    T &make(Args... args)
        { return *newptr<T>(args...); }


    /** Move to allocated location from an rvalue reference */
    template<class T>
    T &move(T &&val)
    {
        T *ptr(new T(std::move(val)));
        deleters.push_back(std::bind(&TmpAlloc::del<T>, ptr));
        return *ptr;
    }


    /** Move to allocated location from an rvalue reference.
    This version is easeier to call than the move() above because
    template parameters can usually be deduced. */
    template<class T>
    void move(T *&ptr, T &&val)
    {
        ptr = new T(std::move(val));
        deleters.push_back(std::bind(&TmpAlloc::del<T>, ptr));
    }


    /** Take over memory management from a unique_ptr.
    Works for all types, no object constructors are required or called.
    @return Raw pointer to the object, now managed by this TmpAlloc. */
    template<class T>
    T *take(std::unique_ptr<T> &&uptr)
    {
        T *ptr(uptr.release());
        deleters.push_back(std::bind(&TmpAlloc::del<T>, ptr));
        return ptr;
    }



    /** Allocate and copy from an existing object */
    template<class T>
    T &copy(T const &val)
        { return *newptr<T>(val); }

    ~TmpAlloc() { free(); }
    void free() {
        for (auto ii(deleters.begin()); ii != deleters.end(); ++ii) {
            (*ii)();
        }
        deleters.clear();
    }

};



}   // namespace ibmisc





#if 0
template <typename T_SRC, typename T_DEST, typename T_DELETER>
bool dynamic_pointer_move(std::unique_ptr<T_DEST, T_DELETER> & dest,
                          std::unique_ptr<T_SRC, T_DELETER> & src)
{
    if (!src) {
        dest.reset();
        return true;
    }

    T_DEST * dest_ptr = dynamic_cast<T_DEST *>(src.get());
    if (!dest_ptr)
        return false;

    std::unique_ptr<T_DEST, T_DELETER> dest_temp(dest_ptr, src.get_deleter());
    src.release();
    dest.swap(dest_temp);
    return true;
}

template <typename T_SRC, typename T_DEST>
bool dynamic_pointer_move(std::unique_ptr<T_DEST> & dest,
                          std::unique_ptr<T_SRC> & src)
{
    if (!src) {
        dest.reset();
        return true;
    }

    T_DEST * dest_ptr = dynamic_cast<T_DEST *>(src.get());
    if (!dest_ptr)
        return false;

    src.release();
    dest.reset(dest_ptr);
    return true;
}
#endif

