use std::cell::RefCell;
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::rc::Rc;
use crate::refs::Pool;

#[derive(Debug, PartialEq, Eq)]
pub struct LightRef<T> {
    inner: Rc<RefCell<T>>,
}
impl<T> Deref for LightRef<T> {
    type Target = Rc<RefCell<T>>;
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}
impl<T> DerefMut for LightRef<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}
impl<T> Clone for LightRef<T> {
    fn clone(&self) -> Self {
        LightRef { inner: self.inner.clone() }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct AllocPool<T> {
    // Not really a pool; instead, track allocations with reference counting
    phantom_data: PhantomData<T>,
}

impl<T> AllocPool<T> {
    pub fn new() -> Self {
        Self {
            phantom_data: PhantomData,
        }
    }
}

impl<T> Pool<T> for AllocPool<T> {
    fn alloc(&self, t: T) -> LightRef<T> {
        LightRef { inner: Rc::new(RefCell::new(t)) }
    }

    fn free(&self, r: LightRef<T>) {
        // Check that the only remaining strong reference to the data is `r`
        assert_eq!(Rc::strong_count(&r.inner), 1, "attempted to free `LightRef` while it still has referents!");

        // When `r` falls out of scope, its reference count will go to 0, and it will be dropped.
    }
}

/// A fake pool that always allocates with `Rc<RefCell<>>`
/// (this is the default behavior when the `fast_unsafe_light_ref` feature is disabled,
/// but is useful to have for tests even when that feature is enabled)
pub type TestPool<T> = AllocPool<T>;
