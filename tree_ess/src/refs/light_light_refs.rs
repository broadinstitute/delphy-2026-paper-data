use std::cell::RefCell;
use std::marker::PhantomData;
use std::ptr::NonNull;
use crate::refs::Pool;

#[derive(Debug)]
pub struct LightRef<T> {
    inner: NonNull<T>,
}
impl<T> LightRef<T> {
    pub fn borrow(&self) -> &T {
        // SAFETY: Caller is responsible to enforce Rust aliasing rules
        // (i.e., at any point in time, either no borrows, 1+ immutable borrows, or 1 mutable borrow)
        unsafe { self.inner.as_ref() }
    }

    pub fn borrow_mut(&self) -> &mut T {
        // SAFETY: Caller is responsible to enforce Rust aliasing rules
        // (i.e., at any point in time, either no borrows, 1+ immutable borrows, or 1 mutable borrow)
        unsafe { NonNull::new_unchecked(self.inner.as_ptr()).as_mut() }
    }
}
impl<T> Clone for LightRef<T> {
    fn clone(&self) -> Self {
        LightRef { inner: self.inner.clone() }
    }
}
impl<T: PartialEq> PartialEq for LightRef<T> {
    fn eq(&self, other: &LightRef<T>) -> bool {
        *self.borrow() == *other.borrow()
    }
}
impl<T: Eq> Eq for LightRef<T> {}

// Objects are allocated on the heap, but memory management is done manually
pub struct AllocPool<T> {
    phantom_data: PhantomData<T>,
}
impl<T> AllocPool<T> {
    pub fn new() -> Self {
        Self { phantom_data: PhantomData }
    }
}
impl<T> Pool<T> for AllocPool<T> {
    // SAFETY: Callers must ensure that no references to objects in the pool remain
    // when the pool is dropped
    fn alloc(&self, t: T) -> LightRef<T> {
        LightRef { inner: unsafe { NonNull::new_unchecked(Box::into_raw(Box::new(t))) } }
    }

    // SAFETY: Manual memory management => you should free allocated objects exactly once
    fn free(&self, r: LightRef<T>) {
        drop(unsafe { Box::from_raw(r.inner.as_ptr()) });
    }
}

/// A fake pool that always allocates with `Rc<RefCell<>>`
/// (this is the default behavior when the `fast_unsafe_light_ref` feature is disabled,
/// but is useful to have for tests even when that feature is enabled)
pub struct TestPool<T> {
    objects: RefCell<Vec<Box<T>>>,
}
impl<T> TestPool<T> {
    pub fn new() -> Self {
        Self { objects: RefCell::new(vec![]) }
    }
}
impl<T> Pool<T> for TestPool<T> {
    // SAFETY: Callers must ensure that no references to objects in the pool remain
    // when the pool is dropped
    fn alloc(&self, t: T) -> LightRef<T> {
        let obj = Box::new(t);
        let ptr = Box::into_raw(obj);
        self.objects.borrow_mut().push(unsafe { Box::from_raw(ptr) });
        LightRef { inner: unsafe { NonNull::new_unchecked(ptr) } }
    }

    fn free(&self, _r: LightRef<T>) {
        // No-op
    }
}

// TODO: Implement other pools, like `FixedSizeArenaPool` or `GrowingPool`