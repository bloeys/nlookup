package nset

import (
	"bytes"
	"fmt"
	"math/bits"
	"reflect"
	"strings"
	"unsafe"
)

var _ fmt.Stringer = &NSet[uint8]{}

type BucketType uint8
type StorageType uint64

const (
	BucketCount = 128
	// StorageTypeBits is the number of bits used per storage unit in each bucket.
	//
	// NOTE: this must be a power of 2, otherwise FastModPower2 will break and must be replaced by a normal x%y
	// NOTE: GetStorageUnitIndex must be adjusted if this value is changed
	StorageTypeBits    = 64
	BucketIndexingBits = 7
)

// IntsIf is limited to uint32 because we can store ALL 4 Billion uint32 numbers
// in 512MB with NSet (instead of the normal 16GB for an array of all uint32s).
// But if we allow uint64 (or int, since int can be 64-bit) users can easily put a big 64-bit number and use more RAM than maybe Google and crash.
type IntsIf interface {
	uint8 | uint16 | uint32
}

type Bucket struct {
	Data             []StorageType
	StorageUnitCount uint32
}

type NSet[T IntsIf] struct {
	Buckets [BucketCount]Bucket
	//StorageUnitCount the number of uint64 integers that are used to indicate presence of numbers in the set
	StorageUnitCount uint32
	shiftAmount      T
	SetBits          uint64
}

func (n *NSet[T]) Add(x T) {

	bucket := n.GetBucketFromValue(x)

	unitIndex := n.GetStorageUnitIndex(x)
	if unitIndex >= bucket.StorageUnitCount {

		storageUnitsToAdd := unitIndex - bucket.StorageUnitCount + 1
		bucket.Data = append(bucket.Data, make([]StorageType, storageUnitsToAdd)...)

		n.StorageUnitCount += storageUnitsToAdd
		bucket.StorageUnitCount += storageUnitsToAdd
	}

	oldStorage := bucket.Data[unitIndex]
	newStorage := oldStorage | n.GetBitMask(x)

	bucket.Data[unitIndex] = newStorage
	n.SetBits += uint64(bits.OnesCount64(uint64(^oldStorage) & uint64(newStorage)))
}

func (n *NSet[T]) AddMany(values ...T) {

	for i := 0; i < len(values); i++ {

		x := values[i]
		bucket := n.GetBucketFromValue(x)

		unitIndex := n.GetStorageUnitIndex(x)
		if unitIndex >= bucket.StorageUnitCount {

			storageUnitsToAdd := unitIndex - bucket.StorageUnitCount + 1
			bucket.Data = append(bucket.Data, make([]StorageType, storageUnitsToAdd)...)

			n.StorageUnitCount += storageUnitsToAdd
			bucket.StorageUnitCount += storageUnitsToAdd
		}

		oldStorage := bucket.Data[unitIndex]
		newStorage := oldStorage | n.GetBitMask(x)

		bucket.Data[unitIndex] = newStorage
		n.SetBits += uint64(bits.OnesCount64(uint64(^oldStorage) & uint64(newStorage)))
	}
}

func (n *NSet[T]) Remove(x T) {

	b := n.GetBucketFromValue(x)
	unitIndex := n.GetStorageUnitIndex(x)
	if unitIndex >= b.StorageUnitCount {
		return
	}

	oldStorage := b.Data[unitIndex]
	newStorage := oldStorage &^ n.GetBitMask(x)

	b.Data[unitIndex] = newStorage
	n.SetBits -= uint64(bits.OnesCount64(uint64(oldStorage) & uint64(^newStorage)))
}

func (n *NSet[T]) Contains(x T) bool {
	return n.isSet(x)
}

func (n *NSet[T]) ContainsAny(values ...T) bool {

	for _, x := range values {
		if n.isSet(x) {
			return true
		}
	}

	return false
}

func (n *NSet[T]) ContainsAll(values ...T) bool {

	for _, x := range values {
		if !n.isSet(x) {
			return false
		}
	}

	return true
}

func (n *NSet[T]) isSet(x T) bool {
	b := n.GetBucketFromValue(x)
	unitIndex := n.GetStorageUnitIndex(x)
	return unitIndex < b.StorageUnitCount && b.Data[unitIndex]&n.GetBitMask(x) != 0
}

func (n *NSet[T]) GetBucketFromValue(x T) *Bucket {
	return &n.Buckets[n.GetBucketIndex(x)]
}

func (n *NSet[T]) GetBucketIndex(x T) BucketType {
	//Use the top 'n' bits as the index to the bucket
	return BucketType(x >> n.shiftAmount)
}

func (n *NSet[T]) GetStorageUnitIndex(x T) uint32 {

	//The top 'n' bits are used to select the bucket so we need to remove them before finding storage
	//unit and bit mask. This is done by shifting left by 4 which removes the top 'n' bits,
	//then shifting right by 4 which puts the bits back to their original place, but now
	//the top 'n' bits are zeros.

	// Since StorageTypeBits is known and is a power of 2, we can replace the division
	// with a right shift.
	//
	// The below return is equal to: return uint32(((x << BucketIndexingBits) >> BucketIndexingBits) / StorageTypeBits)
	return uint32(((x << BucketIndexingBits) >> BucketIndexingBits) >> 6)
}

func (n *NSet[T]) GetBitMask(x T) StorageType {
	//Removes top 'n' bits
	return 1 << FastModPower2(((x<<BucketIndexingBits)>>BucketIndexingBits), StorageTypeBits)
}

func (n *NSet[T]) Union(otherSet *NSet[T]) {

	for i := 0; i < BucketCount; i++ {

		b1 := &n.Buckets[i]
		b2 := &otherSet.Buckets[i]

		if b1.StorageUnitCount < b2.StorageUnitCount {

			storageUnitsToAdd := b2.StorageUnitCount - b1.StorageUnitCount
			b1.Data = append(b1.Data, make([]StorageType, storageUnitsToAdd)...)

			b1.StorageUnitCount += storageUnitsToAdd
			n.StorageUnitCount += storageUnitsToAdd
		}

		for j := 0; j < len(b1.Data) && j < len(b2.Data); j++ {

			oldStorage := b1.Data[j]
			newStorage := oldStorage | b2.Data[j]

			b1.Data[j] = newStorage
			n.SetBits += uint64(bits.OnesCount64(uint64(^oldStorage) & uint64(newStorage)))
		}
	}
}

func (n *NSet[T]) GetIntersection(otherSet *NSet[T]) *NSet[T] {

	outSet := NewNSet[T]()

	for i := 0; i < BucketCount; i++ {

		b1 := &n.Buckets[i]
		b2 := &otherSet.Buckets[i]

		newB := &outSet.Buckets[i]
		for j := uint32(0); j < b1.StorageUnitCount && j < b2.StorageUnitCount; j++ {

			if b1.Data[j]&b2.Data[j] == 0 {
				continue
			}

			if newB.StorageUnitCount < j+1 {
				storageUnitsToAdd := j + 1 - newB.StorageUnitCount
				newB.Data = append(newB.Data, make([]StorageType, storageUnitsToAdd)...)

				newB.StorageUnitCount += storageUnitsToAdd
				outSet.StorageUnitCount += storageUnitsToAdd
			}

			newStorage := b1.Data[j] & b2.Data[j]
			newB.Data[j] = newStorage
			outSet.SetBits += uint64(bits.OnesCount64(uint64(newStorage)))
		}
	}

	return outSet
}

// GetAllElements returns all the added numbers added to NSet.
//
// NOTE: Be careful with this if you have a lot of elements in NSet because NSet is compressed while the returned array is not.
// In the worst case (all uint32s stored) the returned array will be ~4.2 billion elements and will use 16+ GBs of RAM.
func (n *NSet[T]) GetAllElements() []T {

	elements := make([]T, 0, n.SetBits)

	if n.SetBits == 0 {
		return elements
	}

	for i := 0; i < BucketCount; i++ {

		//bucketIndexBits are the bits removed from the original value to use for bucket indexing.
		//We will use this to restore the original value 'x' once an intersection is detected
		bucketIndexBits := T(i << n.shiftAmount)

		b1 := &n.Buckets[i]
		for j := 0; j < len(b1.Data); j++ {

			storageUnit := b1.Data[j]
			if storageUnit == 0 {
				continue
			}

			onesCount := bits.OnesCount64(uint64(storageUnit))

			mask := StorageType(1 << 0)                                     //This will be used to check set bits. Numbers will be reconstructed only for set bits
			firstStorageUnitValue := T(j*StorageTypeBits) | bucketIndexBits //StorageUnitIndex = noBucketBitsX / StorageTypeBits. So: noBucketBitsX = StorageUnitIndex * StorageTypeBits; Then: x = noBucketBitsX | bucketIndexBits

			for k := T(0); onesCount > 0 && k < StorageTypeBits; k++ {

				if storageUnit&mask > 0 {
					elements = append(elements, firstStorageUnitValue+k)
					onesCount--
				}

				mask <<= 1
			}
		}
	}

	return elements
}

func (n *NSet[T]) IsEq(otherSet *NSet[T]) bool {

	if n.SetBits != otherSet.SetBits {
		return false
	}

	//Equal storage unit count doesn't mean all buckets have same size, so we check per bucket
	for i := 0; i < len(n.Buckets); i++ {
		if n.Buckets[i].StorageUnitCount != otherSet.Buckets[i].StorageUnitCount {
			return false
		}
	}

	for i := 0; i < len(n.Buckets); i++ {

		b1 := &n.Buckets[i]
		b2 := &otherSet.Buckets[i]

		bucketsEqual := (b1.StorageUnitCount == 0 && b2.StorageUnitCount == 0) || bytes.Equal(
			unsafe.Slice((*byte)(unsafe.Pointer(&b1.Data[0])), len(b1.Data)*int(unsafe.Sizeof(b1.Data[0]))),
			unsafe.Slice((*byte)(unsafe.Pointer(&b2.Data[0])), len(b2.Data)*int(unsafe.Sizeof(b2.Data[0]))),
		)

		if !bucketsEqual {
			return false
		}
	}

	return true
}

func (n *NSet[T]) HasIntersection(otherSet *NSet[T]) bool {

	for i := 0; i < len(n.Buckets); i++ {

		b1 := &n.Buckets[i]
		b2 := &otherSet.Buckets[i]

		for j := 0; j < len(b1.Data) && j < len(b2.Data); j++ {

			if b1.Data[j]&b2.Data[j] > 0 {
				return true
			}
		}
	}

	return false
}

// String returns a string of the storage as bytes separated by spaces. A comma is between each storage unit
func (n *NSet[T]) String() string {

	b := strings.Builder{}
	b.Grow(int(n.StorageUnitCount*StorageTypeBits + n.StorageUnitCount*2))

	for i := 0; i < len(n.Buckets); i++ {

		bucket := &n.Buckets[i]
		for j := 0; j < len(bucket.Data); j++ {

			x := bucket.Data[j]
			shiftAmount := StorageTypeBits - 8
			for shiftAmount >= 0 {

				byteToShow := uint8(x >> shiftAmount)
				if shiftAmount > 0 {
					b.WriteString(fmt.Sprintf("%08b ", byteToShow))
				} else {
					b.WriteString(fmt.Sprintf("%08b", byteToShow))
				}

				shiftAmount -= 8
			}
			b.WriteString(", ")
		}
	}

	return b.String()
}

func (n *NSet[T]) Copy() *NSet[T] {

	newSet := NewNSet[T]()
	for i := 0; i < len(n.Buckets); i++ {

		b := &n.Buckets[i]
		newB := &newSet.Buckets[i]

		newB.StorageUnitCount = b.StorageUnitCount
		newB.Data = make([]StorageType, len(b.Data))

		copy(newB.Data, b.Data)
	}

	newSet.StorageUnitCount = n.StorageUnitCount
	return newSet

}

// Len returns the number of values stored (i.e. bits set to 1).
// It is the same as NSet.SetBits.
func (n *NSet[T]) Len() uint64 {
	return n.SetBits
}

func UnionSets[T IntsIf](set1, set2 *NSet[T]) *NSet[T] {

	newSet := NewNSet[T]()

	// This is an optimization that makes it so that we only need to count bits
	// when doing union with set2
	newSet.SetBits = set1.SetBits

	for i := 0; i < BucketCount; i++ {

		b1 := &set1.Buckets[i]
		b2 := &set2.Buckets[i]

		//Size bucket
		bucketSize := b1.StorageUnitCount
		if b2.StorageUnitCount > bucketSize {
			bucketSize = b2.StorageUnitCount
		}

		newB := &newSet.Buckets[i]
		newB.Data = make([]StorageType, bucketSize)

		newB.StorageUnitCount = bucketSize
		newSet.StorageUnitCount += bucketSize

		//Union fields of both sets on the new set
		copy(newB.Data, b1.Data)

		for j := 0; j < len(b2.Data); j++ {

			oldStorage := newB.Data[j]
			newStorage := oldStorage | b2.Data[j]

			newB.Data[j] = newStorage
			newSet.SetBits += uint64(bits.OnesCount64(uint64(^oldStorage) & uint64(newStorage)))
		}
	}

	return newSet
}

// FastModPower2 is a fast version of x%y that only works when y is a power of 2
func FastModPower2[T uint8 | uint16 | uint32 | uint64](x, y T) T {
	return x & (y - 1)
}

func NewNSet[T IntsIf]() *NSet[T] {

	n := &NSet[T]{
		Buckets:          [BucketCount]Bucket{},
		StorageUnitCount: 0,
		//We use this to either extract or clear the top 'n' bits, as they are used to select the bucket
		shiftAmount: T(reflect.TypeOf(*new(T)).Bits()) - BucketIndexingBits,
	}

	for i := 0; i < len(n.Buckets); i++ {
		n.Buckets[i].Data = make([]StorageType, 0)
	}

	return n
}
