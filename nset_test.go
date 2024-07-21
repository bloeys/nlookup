package nset_test

import (
	"fmt"
	"math"
	"math/rand"
	"testing"

	"github.com/bloeys/nset"
)

const (
	maxBenchSize = 10_000_000
	RandSeed     = 9_812_938_704
)

var (
	dump          int
	fullRangeNSet *nset.NSet[uint32]
)

func TestNSet(t *testing.T) {

	n1 := nset.NewNSet[uint32]()

	// Double add/remove of the same value is not only important to test SetBits, but also
	// to test for bugs where double adding/removing incorrectly flips bits (checked using Contains())
	n1.Add(0)
	AllTrue(t, n1.Len() == 1)

	n1.Add(0)
	AllTrue(t, n1.Len() == 1)

	n1.Add(1)
	AllTrue(t, n1.Len() == 2)

	n1.Add(63)
	AllTrue(t, n1.Len() == 3)

	n1.Add(math.MaxUint32)
	AllTrue(t, n1.Len() == 4)

	AllTrue(t, n1.Contains(0), n1.Contains(1), n1.Contains(63), n1.Contains(math.MaxUint32), !n1.Contains(10), !n1.Contains(599))
	AllTrue(t, n1.ContainsAll(0, 1, 63), !n1.ContainsAll(9, 0, 1), !n1.ContainsAll(0, 1, 63, 99))
	AllTrue(t, n1.ContainsAny(0, 1, 63), n1.ContainsAny(9, 99, 999, 1), !n1.ContainsAny(9, 99, 999))

	IsEq(t, nset.BucketCount-1, n1.GetBucketIndex(math.MaxUint32))
	IsEq(t, math.MaxUint32/64/nset.BucketCount, n1.GetStorageUnitIndex(math.MaxUint32))

	nCopy := n1.Copy()

	n1.Remove(1)
	AllTrue(t, n1.Len() == 3)

	n1.Remove(1)
	AllTrue(t, n1.Len() == 3)

	AllTrue(t, n1.Contains(0), n1.Contains(63), !n1.Contains(1), nCopy.ContainsAll(0, 1, 63, math.MaxUint32))

	//Intersections
	n2 := nset.NewNSet[uint32]()
	n2.AddMany(1000, 63, 5, 10)

	n3 := nset.NewNSet[uint32]()
	n3.AddMany(math.MaxUint32)

	AllTrue(t, n1.HasIntersection(n2), n2.HasIntersection(n1), n3.HasIntersection(n1), !n3.HasIntersection(n2))

	n4 := nset.NewNSet[uint32]()
	n4.AddMany(0, 1, 64, math.MaxUint32)

	n5 := nset.NewNSet[uint32]()
	n5.AddMany(0, 1, 63, 64, math.MaxUint32)

	n4n5 := n4.GetIntersection(n5)

	n4n5Twin := nset.NewNSet[uint32]()
	n4n5Twin.AddMany(0, 1, 64, math.MaxUint32)

	AllTrue(t, n4n5.Len() == 4, n4n5.Len() == n4n5Twin.Len(), n4n5.ContainsAll(0, 1, 64, math.MaxUint32), !n4n5.Contains(63), n4n5Twin.IsEq(n4n5))

	//Union
	n6 := nset.NewNSet[uint32]()
	n6.AddMany(1, 4, 7, 100, 1000)

	n7 := nset.NewNSet[uint32]()
	n7.AddMany(1, math.MaxUint32)
	n7OldStorageUnitCount := n7.StorageUnitCount
	n7.Union(n6)

	AllTrue(t, n6.Len() == 5, n7.Len() == 6, n6.ContainsAll(1, 4, 7, 100, 1000), !n6.Contains(math.MaxUint32), n7.ContainsAll(1, 4, 7, 100, 1000, math.MaxUint32), n7.StorageUnitCount == n7OldStorageUnitCount+n6.StorageUnitCount-1)

	//UnionSets
	n7 = nset.NewNSet[uint32]()
	n7.AddMany(4, math.MaxUint32)

	unionedSet := nset.UnionSets(n6, n7)
	AllTrue(t, unionedSet.Len() == 6, !n6.Contains(math.MaxUint32), !n7.ContainsAny(7, 100, 1000), unionedSet.ContainsAll(4, 7, 100, 1000, math.MaxUint32), unionedSet.StorageUnitCount == n6.StorageUnitCount+n7OldStorageUnitCount-1)

	//Equality
	AllTrue(t, !n6.IsEq(n7))

	n7.Union(n6)
	AllTrue(t, !n6.IsEq(n7))

	n6.Union(n7)
	AllTrue(t, n6.IsEq(n7))

	//GetAllElements
	n8 := nset.NewNSet[uint32]()
	n8.AddMany(0, 1, 55, 1000, 10000)

	n8Elements := n8.GetAllElements()
	AllTrue(t, len(n8Elements) == 5, n8Elements[0] == 0, n8Elements[1] == 1, n8Elements[2] == 55, n8Elements[3] == 1000, n8Elements[4] == 10000)
}

func TestNSetFullRange(t *testing.T) {

	if fullRangeNSet == nil {

		fullRangeNSet = nset.NewNSet[uint32]()
		println("Adding all uint32 to NSet...")
		for i := uint32(0); i < math.MaxUint32; i++ {
			fullRangeNSet.Add(i)
			if i%1_000_000_000 == 0 {
				fmt.Printf("i=%d billion\n", i/1_000_000_000)
			}
		}
		fullRangeNSet.Add(math.MaxUint32)
		AllTrue(t, fullRangeNSet.SetBits == math.MaxUint32+1)
	}

	n := fullRangeNSet
	IsEq(t, 67_108_864, n.StorageUnitCount)
	for i := 0; i < len(n.Buckets); i++ {

		b := &n.Buckets[i]
		IsEq(t, 524288, b.StorageUnitCount)

		for j := 0; j < len(b.Data); j++ {
			if b.Data[j] != math.MaxUint64 {
				t.Fatalf("Error: storage unit is NOT equal to MaxUint64 (i=%d,j=%d)! Expected math.MaxUint64 but got '%08b'\n",
					i,
					j,
					b.Data[j])
			}
		}
	}

}

func AllTrue(t *testing.T, values ...bool) (success bool) {

	success = true

	for i := 0; i < len(values); i++ {
		if !values[i] {
			t.Fatalf("Expected 'true' but got 'false'\n")
			success = false
		}
	}

	return success
}

func IsEq[T comparable](t *testing.T, expected, val T) bool {

	if val == expected {
		return true
	}

	t.Fatalf("Expected '%v' but got '%v'\n", expected, val)
	return false
}

func BenchmarkNSetAdd(b *testing.B) {

	n := nset.NewNSet[uint32]()

	for i := uint32(0); i < uint32(b.N); i++ {
		n.Add(i % maxBenchSize)
	}
}

func BenchmarkMapAdd(b *testing.B) {

	hMap := map[uint32]struct{}{}

	for i := uint32(0); i < uint32(b.N); i++ {
		hMap[i%maxBenchSize] = struct{}{}
	}
}

func BenchmarkNSetAddRand(b *testing.B) {

	n := nset.NewNSet[uint32]()

	rand := rand.New(rand.NewSource(RandSeed))
	for i := 0; i < b.N; i++ {
		n.Add(rand.Uint32() % maxBenchSize)
	}
}

func BenchmarkNSetAddRandNoSizeLimit(b *testing.B) {

	n := nset.NewNSet[uint32]()

	rand := rand.New(rand.NewSource(RandSeed))
	for i := 0; i < b.N; i++ {
		n.Add(rand.Uint32())
	}
}

func BenchmarkMapAddRand(b *testing.B) {

	hMap := map[uint32]struct{}{}

	rand := rand.New(rand.NewSource(RandSeed))
	for i := 0; i < b.N; i++ {
		hMap[rand.Uint32()%maxBenchSize] = struct{}{}
	}
}

func BenchmarkNSetContains(b *testing.B) {

	//Init
	b.StopTimer()
	n := nset.NewNSet[uint32]()

	for i := uint32(0); i < maxBenchSize; i++ {
		n.Add(i)
	}
	b.StartTimer()

	//Work
	found := 0
	for i := uint32(0); i < uint32(b.N); i++ {
		if n.Contains(i) {
			found++
		}
	}

	dump = found
}

func BenchmarkMapContains(b *testing.B) {

	//Init
	b.StopTimer()
	hMap := map[uint32]struct{}{}

	for i := uint32(0); i < maxBenchSize; i++ {
		hMap[i] = struct{}{}
	}
	b.StartTimer()

	//Work
	found := 0
	for i := uint32(0); i < uint32(b.N); i++ {
		if _, ok := hMap[i]; ok {
			found++
		}
	}

	dump = found
}

func BenchmarkNSetContainsRand(b *testing.B) {

	//Init
	b.StopTimer()
	n := nset.NewNSet[uint32]()

	for i := uint32(0); i < maxBenchSize; i++ {
		n.Add(i)
	}
	b.StartTimer()

	//Work
	found := 0
	rand := rand.New(rand.NewSource(RandSeed))
	for i := 0; i < b.N; i++ {

		randVal := rand.Uint32()
		if n.Contains(randVal) {
			found++
		}
	}

	dump = found
}

func BenchmarkNSetContainsRandFullRange(b *testing.B) {

	//Init
	if fullRangeNSet == nil {

		b.StopTimer()

		fullRangeNSet = nset.NewNSet[uint32]()
		println("Preparing full range NSet...")
		for i := uint32(0); i < math.MaxUint32; i++ {
			fullRangeNSet.Add(i)
		}
		fullRangeNSet.Add(math.MaxUint32)

		b.StartTimer()
	}

	n := fullRangeNSet

	//Work
	found := 0
	rand := rand.New(rand.NewSource(RandSeed))
	for i := 0; i < b.N; i++ {

		randVal := rand.Uint32()
		if n.Contains(randVal) {
			found++
		}
	}

	dump = found
}

func BenchmarkMapContainsRand(b *testing.B) {

	//Init
	b.StopTimer()
	hMap := map[uint32]struct{}{}

	for i := uint32(0); i < maxBenchSize; i++ {
		hMap[i] = struct{}{}
	}
	b.StartTimer()

	//Work
	found := 0
	rand := rand.New(rand.NewSource(RandSeed))
	for i := 0; i < b.N; i++ {

		randVal := rand.Uint32()
		if _, ok := hMap[randVal]; ok {
			found++
		}
	}

	dump = found
}

func BenchmarkNSetDelete(b *testing.B) {

	//Init
	b.StopTimer()
	n := nset.NewNSet[uint32]()

	for i := uint32(0); i < maxBenchSize; i++ {
		n.Add(i)
	}
	b.StartTimer()

	//Work
	for i := uint32(0); i < uint32(b.N); i++ {
		n.Remove(i)
	}
}

func BenchmarkMapDelete(b *testing.B) {

	//Init
	b.StopTimer()
	hMap := map[uint32]struct{}{}

	for i := uint32(0); i < maxBenchSize; i++ {
		hMap[i] = struct{}{}
	}
	b.StartTimer()

	//Work
	for i := uint32(0); i < uint32(b.N); i++ {
		delete(hMap, i)
	}
}

func BenchmarkNSetDeleteRand(b *testing.B) {

	//Init
	b.StopTimer()
	n := nset.NewNSet[uint32]()

	for i := uint32(0); i < maxBenchSize; i++ {
		n.Add(i)
	}
	b.StartTimer()

	//Work
	rand := rand.New(rand.NewSource(RandSeed))
	for i := 0; i < b.N; i++ {

		randVal := rand.Uint32()
		n.Remove(randVal)
	}
}

func BenchmarkMapDeleteRand(b *testing.B) {

	//Init
	b.StopTimer()
	hMap := map[uint32]struct{}{}

	for i := uint32(0); i < maxBenchSize; i++ {
		hMap[i] = struct{}{}
	}
	b.StartTimer()

	//Work
	rand := rand.New(rand.NewSource(RandSeed))
	for i := 0; i < b.N; i++ {

		randVal := rand.Uint32()
		delete(hMap, randVal)
	}
}

func BenchmarkNSetIsEq(b *testing.B) {

	b.StopTimer()
	s1 := nset.NewNSet[uint32]()
	s2 := nset.NewNSet[uint32]()
	for i := uint32(0); i < maxBenchSize; i++ {
		s1.Add(i)
		s2.Add(i)
	}
	b.StartTimer()

	for i := 0; i < b.N; i++ {
		s1.IsEq(s2)
	}
}

func BenchmarkMapIsEq(b *testing.B) {

	b.StopTimer()
	m1 := map[uint32]struct{}{}
	m2 := map[uint32]struct{}{}
	for i := uint32(0); i < maxBenchSize; i++ {
		m1[i] = struct{}{}
		m2[i] = struct{}{}
	}
	b.StartTimer()

	mapsAreEq := func(m1, m2 map[uint32]struct{}) bool {

		if len(m1) != len(m2) {
			return false
		}

		for k := range m1 {
			if _, ok := m2[k]; !ok {
				return false
			}
		}

		return true
	}

	for i := 0; i < b.N; i++ {
		mapsAreEq(m1, m2)
	}
}

func BenchmarkNSetIsEqRand(b *testing.B) {

	b.StopTimer()

	rand := rand.New(rand.NewSource(RandSeed))
	s1 := nset.NewNSet[uint32]()
	s2 := nset.NewNSet[uint32]()
	for i := uint32(0); i < maxBenchSize; i++ {
		r := rand.Uint32()
		s1.Add(r)
		s2.Add(r)
	}
	b.StartTimer()

	for i := 0; i < b.N; i++ {
		s1.IsEq(s2)
	}
}

func BenchmarkMapIsEqRand(b *testing.B) {

	b.StopTimer()

	rand := rand.New(rand.NewSource(RandSeed))
	m1 := map[uint32]struct{}{}
	m2 := map[uint32]struct{}{}
	for i := uint32(0); i < maxBenchSize; i++ {
		r := rand.Uint32()
		m1[r] = struct{}{}
		m2[r] = struct{}{}
	}
	b.StartTimer()

	mapsAreEq := func(m1, m2 map[uint32]struct{}) bool {

		if len(m1) != len(m2) {
			return false
		}

		for k := range m1 {
			if _, ok := m2[k]; !ok {
				return false
			}
		}

		return true
	}

	for i := 0; i < b.N; i++ {
		mapsAreEq(m1, m2)
	}
}

func BenchmarkNSetIsEqRand100Mil(b *testing.B) {

	b.StopTimer()

	rand := rand.New(rand.NewSource(RandSeed))
	s1 := nset.NewNSet[uint32]()
	s2 := nset.NewNSet[uint32]()
	for i := uint32(0); i < 100_000_000; i++ {
		r := rand.Uint32()
		s1.Add(r)
		s2.Add(r)
	}
	b.StartTimer()

	for i := 0; i < b.N; i++ {
		s1.IsEq(s2)
	}
}

func BenchmarkMapIsEqRand100Mil(b *testing.B) {

	b.StopTimer()

	rand := rand.New(rand.NewSource(RandSeed))
	m1 := map[uint32]struct{}{}
	m2 := map[uint32]struct{}{}
	for i := uint32(0); i < 100_000_000; i++ {
		r := rand.Uint32()
		m1[r] = struct{}{}
		m2[r] = struct{}{}
	}
	b.StartTimer()

	mapsAreEq := func(m1, m2 map[uint32]struct{}) bool {

		if len(m1) != len(m2) {
			return false
		}

		for k := range m1 {
			if _, ok := m2[k]; !ok {
				return false
			}
		}

		return true
	}

	for i := 0; i < b.N; i++ {
		mapsAreEq(m1, m2)
	}
}

var getIntersectionNset *nset.NSet[uint32]

func BenchmarkNSetGetIntersection(b *testing.B) {

	b.StopTimer()
	s1 := nset.NewNSet[uint32]()
	s2 := nset.NewNSet[uint32]()
	for i := uint32(0); i < maxBenchSize; i++ {
		s1.Add(i)
		s2.Add(i)
	}
	b.StartTimer()

	for i := 0; i < b.N; i++ {
		getIntersectionNset = s1.GetIntersection(s2)
	}
}

var getIntersectionTempMap map[uint32]struct{}

func BenchmarkMapGetIntersection(b *testing.B) {

	b.StopTimer()
	m1 := map[uint32]struct{}{}
	m2 := map[uint32]struct{}{}
	for i := uint32(0); i < maxBenchSize; i++ {
		m1[i] = struct{}{}
		m2[i] = struct{}{}
	}
	b.StartTimer()

	getIntersection := func(m1, m2 map[uint32]struct{}) map[uint32]struct{} {

		outMap := map[uint32]struct{}{}

		for k := range m1 {
			if _, ok := m2[k]; ok {
				outMap[k] = struct{}{}
			}
		}

		return outMap
	}

	for i := 0; i < b.N; i++ {
		getIntersectionTempMap = getIntersection(m1, m2)
	}
}

func BenchmarkNSetGetIntersectionRand(b *testing.B) {

	b.StopTimer()

	rand := rand.New(rand.NewSource(RandSeed))

	s1 := nset.NewNSet[uint32]()
	s2 := nset.NewNSet[uint32]()
	for i := uint32(0); i < maxBenchSize; i++ {

		r := rand.Uint32()
		s1.Add(r)
		s2.Add(r)
	}
	b.StartTimer()

	for i := 0; i < b.N; i++ {
		getIntersectionNset = s1.GetIntersection(s2)
	}
}

func BenchmarkMapGetIntersectionRand(b *testing.B) {

	b.StopTimer()

	rand := rand.New(rand.NewSource(RandSeed))

	m1 := map[uint32]struct{}{}
	m2 := map[uint32]struct{}{}
	for i := uint32(0); i < maxBenchSize; i++ {

		r := rand.Uint32()
		m1[r] = struct{}{}
		m2[r] = struct{}{}
	}
	b.StartTimer()

	getIntersection := func(m1, m2 map[uint32]struct{}) map[uint32]struct{} {

		outMap := map[uint32]struct{}{}

		for k := range m1 {
			if _, ok := m2[k]; ok {
				outMap[k] = struct{}{}
			}
		}

		return outMap
	}

	for i := 0; i < b.N; i++ {
		getIntersectionTempMap = getIntersection(m1, m2)
	}
}

var elementCount int

func BenchmarkNSetGetAllElements(b *testing.B) {

	b.StopTimer()

	s1 := nset.NewNSet[uint32]()
	for i := uint32(0); i < maxBenchSize; i++ {
		s1.Add(i)
	}
	b.StartTimer()

	var elements []uint32
	for i := 0; i < b.N; i++ {
		elements = s1.GetAllElements()
	}

	elementCount = len(elements)
}

func BenchmarkMapGetAllElements(b *testing.B) {

	b.StopTimer()

	m1 := map[uint32]struct{}{}
	for i := uint32(0); i < maxBenchSize; i++ {
		m1[i] = struct{}{}
	}
	b.StartTimer()

	getElementsFunc := func(m map[uint32]struct{}) []uint32 {

		e := make([]uint32, 0, len(m))
		for k := range m {
			e = append(e, k)
		}

		return e
	}

	var elements []uint32
	for i := 0; i < b.N; i++ {
		elements = getElementsFunc(m1)
	}

	elementCount = len(elements)
}

func BenchmarkNSetGetAllElementsRand(b *testing.B) {

	b.StopTimer()

	rand := rand.New(rand.NewSource(RandSeed))
	s1 := nset.NewNSet[uint32]()
	for i := uint32(0); i < maxBenchSize; i++ {
		s1.Add(rand.Uint32())
	}
	b.StartTimer()

	var elements []uint32
	for i := 0; i < b.N; i++ {
		elements = s1.GetAllElements()
	}

	elementCount = len(elements)
}

func BenchmarkMapGetAllElementsRand(b *testing.B) {

	b.StopTimer()

	rand := rand.New(rand.NewSource(RandSeed))

	m1 := map[uint32]struct{}{}
	for i := uint32(0); i < maxBenchSize; i++ {
		m1[rand.Uint32()] = struct{}{}
	}

	getElementsFunc := func(m map[uint32]struct{}) []uint32 {

		e := make([]uint32, 0, len(m))
		for k := range m {
			e = append(e, k)
		}

		return e
	}
	b.StartTimer()

	var elements []uint32
	for i := 0; i < b.N; i++ {
		elements = getElementsFunc(m1)
	}

	elementCount = len(elements)
}

var unionSize int

func BenchmarkNSetUnion(b *testing.B) {

	b.StopTimer()

	s1 := nset.NewNSet[uint32]()
	s2 := nset.NewNSet[uint32]()
	for i := uint32(0); i < maxBenchSize; i++ {
		s1.Add(i)
		s2.Add(i)
	}
	b.StartTimer()

	var union *nset.NSet[uint32]
	for i := 0; i < b.N; i++ {
		union = nset.UnionSets(s1, s2)
	}

	unionSize = int(union.StorageUnitCount)
}

func BenchmarkMapUnion(b *testing.B) {

	b.StopTimer()

	m1 := map[uint32]struct{}{}
	m2 := map[uint32]struct{}{}
	for i := uint32(0); i < maxBenchSize; i++ {
		m1[i] = struct{}{}
		m2[i] = struct{}{}
	}
	b.StartTimer()

	unionFunc := func(m1, m2 map[uint32]struct{}) map[uint32]struct{} {

		u := make(map[uint32]struct{}, len(m1))
		for k := range m1 {
			u[k] = struct{}{}
		}

		for k := range m2 {
			u[k] = struct{}{}
		}

		return u
	}

	var union map[uint32]struct{}
	for i := 0; i < b.N; i++ {
		union = unionFunc(m1, m2)
	}

	unionSize = len(union)
}

func BenchmarkNSetUnionRand(b *testing.B) {

	b.StopTimer()

	randGen := rand.New(rand.NewSource(RandSeed))

	s1 := nset.NewNSet[uint32]()
	s2 := nset.NewNSet[uint32]()
	for i := uint32(0); i < maxBenchSize; i++ {
		r := randGen.Uint32()
		s1.Add(r)
		s2.Add(r)
	}
	b.StartTimer()

	var union *nset.NSet[uint32]
	for i := 0; i < b.N; i++ {
		union = nset.UnionSets(s1, s2)
	}

	unionSize = int(union.StorageUnitCount)
}

func BenchmarkMapUnionRand(b *testing.B) {

	b.StopTimer()

	randGen := rand.New(rand.NewSource(RandSeed))

	m1 := map[uint32]struct{}{}
	m2 := map[uint32]struct{}{}
	for i := uint32(0); i < maxBenchSize; i++ {
		r := randGen.Uint32()
		m1[r] = struct{}{}
		m2[r] = struct{}{}
	}
	b.StartTimer()

	unionFunc := func(m1, m2 map[uint32]struct{}) map[uint32]struct{} {

		u := make(map[uint32]struct{}, len(m1))
		for k := range m1 {
			u[k] = struct{}{}
		}

		for k := range m2 {
			u[k] = struct{}{}
		}

		return u
	}

	var union map[uint32]struct{}
	for i := 0; i < b.N; i++ {
		union = unionFunc(m1, m2)
	}

	unionSize = len(union)
}
