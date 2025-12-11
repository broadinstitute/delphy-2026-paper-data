//! References that behave like `Rc<RefCell<T>>` but where we can trade safety for zero overhead
//! at compile time.
//!
//! The intent is that for highly-linked data structures, once we have verified that all traversals
//! and manipulations respect Rust's aliasing rules, it's very beneficial to get rid of the
//! runtime borrow-checking enforced by [`RefCell<T>`].  Moreover, in such structures, we usually
//! have clear knowledge of when objects are no longer referenced and can be dropped manually,
//! which opens up the possibility of removing the overhead of [`Rc<T>`] as well.
//!
//! [`LightRef<T>`] behaves as a [`Rc<RefCell<T>>`], and is in fact a thin wrapper around it
//! in "debug" mode.  In "production" mode, it is instead a wrapper around [`NonNull<T>`].
//! Since the referent is not reference-counted in "production" mode, allocations and deallocations
//! must be managed manually via an object that implements the [`Pool<T>`] interface.  In "debug"
//! mode, [`AllocPool<T>`] simply creates [`Rc<RefCell<T>>`] values on allocation, and checks that
//! (a) on deallocation, no references to the value remain; and (b) when the pool itself is dropped,
//! all previously allocated objects have been deallocated.  If users of `AllocPool<T>` behave in
//! a way that meets these two criteria, other pool implementations may be used with less overhead.
//! For example, [`FixedSizeArenaPool<T>`] preallocates a fixed-size block of memory from which
//! allocations can be served directly.

#[cfg(not(feature = "fast_unsafe_light_refs"))]
mod heavy_light_refs;

#[cfg(feature = "fast_unsafe_light_refs")]
mod light_light_refs;

#[cfg(not(feature = "fast_unsafe_light_refs"))]
pub use heavy_light_refs::{LightRef, AllocPool, TestPool};

#[cfg(feature = "fast_unsafe_light_refs")]
pub use light_light_refs::{LightRef, AllocPool, TestPool};

pub trait Pool<T> {
    fn alloc(&self, t: T) -> LightRef<T>;
    fn free(&self, r: LightRef<T>);
}