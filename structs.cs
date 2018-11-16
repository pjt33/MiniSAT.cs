/******************************************************************************************
MiniSat -- Copyright (c) 2003-2005, Niklas Een, Niklas Sorensson
MiniSatCS -- Copyright (c) 2006, Michal Moskal

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

using System;
using System.IO;
using System.Text;
using System.Diagnostics;
using System.Collections.Generic;

namespace MiniSatCS {
//=================================================================================================
// 'vec' -- automatically resizable arrays (via 'push()' method):

public class vec<T> {
    T[]  data;
    int sz;

    void     grow(int min_cap) {
		int cap = data == null ? 0 : data.Length;
		if (min_cap <= cap) return;
		if (cap == 0) cap = (min_cap >= 2) ? min_cap : 2;
		else          do cap = (cap*3+1) >> 1; while (cap < min_cap);
		Array.Resize(ref data, cap);
    }

    // Constructors:
    public vec()                { }
    public vec(int size)            { growTo(size); }
    public vec(int size, T pad)     { growTo(size, pad); }
	// (takes ownership of array)
    public vec(T[] array)          
	{
	  data = array;
	  sz = data.Length;
	}

    // Size operations:
    public int      size   ()       { return sz; }
    public void     shrink (int nelems)       
	{
		Solver.assert(nelems <= sz);
		Solver.assert(nelems >= 0);
		sz -= nelems;
	}
	public void     shrinkTo (int nelems)
	{
		shrink (size() - nelems);
	}
    public void     pop    ()             { Solver.assert (sz > 0); sz--; }
    public void     growTo (int size)
	{
		grow(size);
		sz = size;
	}

    public void growTo (int size, T pad)
	{
		int i = sz;

		growTo(size);

		for (; i < sz; ++i)
			data[i] = pad;
	}

    public void     clear  ()
	{
		sz = 0;
	}

    public void     clear  (bool dealloc)
	{
		sz = 0;
	}

    public void     capacity (int size) { grow(size); }

    // Stack interface:
    public void     push  (T elem)     { if (data == null || sz == data.Length) grow(sz+1); data[sz++] = elem; }
    public T       last  ()              { return data[sz-1]; }

    // Vector interface:
	public T this[int index]
	{
		get {
			Solver.assert(index < sz);
			return data[index];
		}
		
		set {
			Solver.assert(index < sz);
			data[index] = value;
		}
	}

    // Duplicatation (preferred instead):
    public void copyTo(vec<T> copy) { copy.clear(); copy.growTo(sz); 
		if (sz > 0) Array.Copy(data, copy.data, sz); }
    public void moveTo(vec<T> dest) { dest.clear(); dest.data = data; dest.sz = sz; 
		data = null; sz = 0;  }
	
	public void sort(IComparer<T> cmp)
	{
		if (sz > 1)
		Array.Sort(data, 0, sz, cmp);
	}

	public IEnumerator<T> GetEnumerator()
	{
		for (int i = 0; i < sz; ++i)
			yield return data[i];
	}

	public void swap(int i, int j)
	{
		Solver.assert(i < sz && j < sz);
		T tmp = data[i];
		data[i] = data[j];
		data[j] = tmp;
	}

	public T[] ToArray()
	{
		T[] res = new T[sz];
		for (int i = 0; i < sz; ++i)
			res[i] = data[i];
		return res;
	}
}


//=================================================================================================


public delegate bool IntLess(int i1, int i2);

public class Heap {
    public IntLess    comp;
    public vec<int> heap = new vec<int>();     // heap of ints
    public vec<int> indices = new vec<int>();  // int . index in heap

	static int left  (int i) { return i+i; }
	static int right (int i) { return i+i + 1; }
	static int parent(int i) { return i >> 1; }

    public void percolateUp(int i)
    {
        int x = heap[i];
        while (parent(i) != 0 && comp(x,heap[parent(i)])){
            heap[i]          = heap[parent(i)];
            indices[heap[i]] = i;
            i                = parent(i);
        }
        heap   [i] = x;
        indices[x] = i;
    }

    public void percolateDown(int i)
    {
        int x = heap[i];
        while (left(i) < heap.size()){
            int child = right(i) < heap.size() && comp(heap[right(i)],heap[left(i)]) ? right(i) : left(i);
            if (!comp(heap[child],x)) break;
            heap[i]          = heap[child];
            indices[heap[i]] = i;
            i                = child;
        }
        heap   [i] = x;
        indices[x] = i;
    }

    public bool ok(int n) { return n >= 0 && n < (int)indices.size(); }

    public Heap(IntLess c) { comp = c; heap.push(-1); }

    public void setBounds (int size) { Solver.assert(size >= 0); indices.growTo(size,0); }
    public bool inHeap    (int n)    { Solver.assert(ok(n)); return indices[n] != 0; }
    public void increase  (int n)    { Solver.assert(ok(n)); Solver.assert(inHeap(n)); percolateUp(indices[n]); }
    public bool empty     ()         { return heap.size() == 1; }

    public void insert(int n) {
		//reportf("{0} {1}\n", n, indices.size());
        Solver.assert(ok(n));
        indices[n] = heap.size();
        heap.push(n);
        percolateUp(indices[n]); }

    public int  getmin() {
        int r            = heap[1];
        heap[1]          = heap.last();
        indices[heap[1]] = 1;
        indices[r]       = 0;
        heap.pop();
        if (heap.size() > 1)
            percolateDown(1);
        return r; }

    public bool heapProperty() {
        return heapProperty(1); }

    public bool heapProperty(int i) {
        return i >= heap.size()
            || ((parent(i) == 0 || !comp(heap[i],heap[parent(i)])) && heapProperty(left(i)) && heapProperty(right(i))); }
}


} // end namespace MiniSat
