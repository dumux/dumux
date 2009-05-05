#ifndef DUNE_ONEDINNDGRID_LIST_HH
#define DUNE_ONEDINNDGRID_LIST_HH

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/iteratorfacades.hh>

namespace Dune {
/** \file
    \brief A simple doubly-linked list needed in OneDInNDGrid
    \todo I'd love to get rid of this and use std::list instead.
    Unfortunately, there are problems.  I need to store pointers/iterators
    within one element which point to another element (e.g. the element father).
 */
    template<class T>
    class OneDInNDGridListIterator 
        : public BidirectionalIteratorFacade<OneDInNDGridListIterator<T>,T>
    {
    public:
        bool equals(const OneDInNDGridListIterator& other) const {
            return pointer_ == other.pointer_;
        }

        T& dereference() {
            return *pointer_;
        }

        void increment() {
            pointer_ = pointer_->succ_;
        }

        void decrement() {
            pointer_ = pointer_->pred_;
        }

        OneDInNDGridListIterator() {}

        OneDInNDGridListIterator(T* pointer) {
            pointer_ = pointer;
        }

        OneDInNDGridListIterator operator=(T* pointer) {
            pointer_ = pointer;
        }

        operator T*() {return pointer_;}

    private:
        T* pointer_;
    };

    template<class T>
    class OneDInNDGridList
    {

    public:

        //typedef OneDInNDGridListIterator<T> iterator;
        typedef T* iterator;
        typedef const T* const_iterator;

        OneDInNDGridList() : numelements(0), begin_(0), rbegin_(0) {}

#if 0
        ~OneDInNDGridList() {
            // Delete all elements
            iterator e = begin();

            while (e) {
            
                iterator eSucc = e->succ_;
                erase(e);
                e = eSucc;
                return;
            }

        }
#endif

        int size() const {return numelements;}
        
        iterator push_back (const T& value) {

            T* i = rbegin();

            // New list element by copy construction
            T* t = new T(value);

            // einfuegen
            if (begin_==0) {
                // einfuegen in leere Liste
                begin_ = t; 
                rbegin_ = t;
            }
            else 
                {
                    // nach Element i.p einsetzen
                    t->pred_ = i;
                    t->succ_ = i->succ_;
                    i->succ_ = t;

                    if (t->succ_!=0) 
                        t->succ_->pred_ = t;

                    // tail neu ?
                    if (rbegin_==i) 
                        rbegin_ = t;
                }

            // Groesse und Rueckgabeiterator
            numelements = numelements+1;

            return t;
        }

        iterator push_back (iterator t) DUNE_DEPRECATED {

            T* i = rbegin();

            // einfuegen
            if (begin_==0) {
                // einfuegen in leere Liste
                begin_ = t; 
                rbegin_ = t;
            }
            else 
                {
                    // nach Element i.p einsetzen
                    t->pred_ = i;
                    t->succ_ = i->succ_;
                    i->succ_ = t;

                    if (t->succ_!=0) 
                        t->succ_->pred_ = t;

                    // tail neu ?
                    if (rbegin_==i) 
                        rbegin_ = t;
                }

            // Groesse und Rueckgabeiterator
            numelements = numelements+1;

            return t;
        }

        iterator insert (iterator i, const T& value) {

            // Insert before 'one-after-the-end' --> append to the list
            if (i==end())
                return push_back(value);

            // New list element by copy construction
            T* t = new T(value);

            // einfuegen
            if (begin_==0) 
                {
                    // einfuegen in leere Liste
                    begin_=t; 
                    rbegin_=t;
                }
            else 
                {
                    // vor Element i.p einsetzen
                    t->succ_ = i;
                    t->pred_ = i->pred_;
                    i->pred_ = t;

                    if (t->pred_!=0) 
                        t->pred_->succ_ = t;
                    // head neu ?
                    if (begin_==i) 
                        begin_ = t;
                }
            
            // Groesse und Rueckgabeiterator
            numelements = numelements+1;
            return t;
        }

        iterator insert (iterator i, iterator t) DUNE_DEPRECATED {

            // Insert before 'one-after-the-end' --> append to the list
            if (i==end())
                return push_back(t);
            
            // einfuegen
            if (begin_==0) 
                {
                    // einfuegen in leere Liste
                    begin_=t; 
                    rbegin_=t;
                }
            else 
                {
                    // vor Element i.p einsetzen
                    t->succ_ = i;
                    t->pred_ = i->pred_;
                    i->pred_ = t;

                    if (t->pred_!=0) 
                        t->pred_->succ_ = t;
                    // head neu ?
                    if (begin_==i) 
                        begin_ = t;
                }
            
            // Groesse und Rueckgabeiterator
            numelements = numelements+1;
            return t;
        }

        void erase (iterator& i)
        {
            // Teste Eingabe
            if (i==0)
                return;
            
            // Ausfaedeln
            if (i->succ_!=0) 
                i->succ_->pred_ = i->pred_;
            if (i->pred_!=0) 
                i->pred_->succ_ = i->succ_;
            
            // head & tail
            if (begin_==i) 
                begin_=i->succ_;
            if (rbegin_==i) 
                rbegin_ = i->pred_;
            
            // Groesse
            numelements = numelements-1;

            // Actually delete the object
            delete(i);
        }

        iterator begin() {
            return begin_;
        }

        const_iterator begin() const {
            return begin_;
        }

        iterator end() {
            return NULL;
        }

        const_iterator end() const {
            return NULL;
        }

        iterator rbegin() {
            return rbegin_;
        }

        const_iterator rbegin() const {
            return rbegin_;
        }

    private:

        int numelements;

        T* begin_;
        T* rbegin_;

    }; // end class OneDInNDGridList
    
} // namespace Dune

#endif
