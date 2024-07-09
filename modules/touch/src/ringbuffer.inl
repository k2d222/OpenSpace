#include <modules/touch/include/ringbuffer.h>

#include <vector>
#include <cstddef> // for size_t

namespace openspace {

template <typename T>
RingBuffer<T>::RingBuffer(size_t capacity)
    : _buffer(capacity)
    , _head(0)
    , _capacity(capacity)
{}

template <typename T>
void RingBuffer<T>::resize(size_t capacity) {
    _buffer.resize(capacity);
    _head = 0;
    _capacity = capacity;
}

template <typename T>
size_t RingBuffer<T>::capacity() const {
    return _capacity;
}

template <typename T>
void RingBuffer<T>::push(const T& value) {
    if (_buffer.size() < _capacity) {
        _buffer.push_back(value);
    } else {
        _buffer[_head] = value;
        _head = (_head + 1) % _capacity;
    }
}

template <typename T>
void RingBuffer<T>::clear() {
    _head = 0;
    _buffer.clear();
}

template <typename T>
size_t RingBuffer<T>::size() const {
    return _buffer.size();
}

template <typename T>
T& RingBuffer<T>::at(size_t index) {
    return _buffer.at((_head + index) % _capacity);
}
template <typename T>
const T& RingBuffer<T>::at(size_t index) const {
    return _buffer.at((_head + index) % _capacity);
}

template <typename T>
T& RingBuffer<T>::front() {
    return _buffer.at(_head);
}
template <typename T>
const T& RingBuffer<T>::front() const {
    return _buffer.at(_head);
}

template <typename T>
T& RingBuffer<T>::back() {
    return _buffer.at((_head + _buffer.size() - 1) % _capacity);
}
template <typename T>
const T& RingBuffer<T>::back() const {
    return _buffer.at((_head + _buffer.size() - 1) % _capacity);
}

} // openspace namespace
