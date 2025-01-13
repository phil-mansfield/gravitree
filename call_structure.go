package gravitree

type CallStructure struct {
	Quantity Quantity
	Sizes1, Sizes2 []int
	NodeDepths []int // For tree1

	// One call for each node in t2.
	TwoSidedLeafCalls [][]Call
	OneSidedLeafCalls [][]Call
	ApproximateCalls [][]Call
}
var _ Quantity = &CallStructure{ } // type-checking

type Call struct {
	I1, I2 int
}

func NewCallStructure(t1, t2 *Tree, q Quantity) *CallStructure {

	sizes1, sizes2 := make([]int, len(t1.Nodes)), make([]int, len(t2.Nodes))
	for i := range t1.Nodes {
		sizes1[i] = t1.Nodes[i].End - t1.Nodes[i].Start
	}
	for i := range t2.Nodes {
		sizes2[i] = t2.Nodes[i].End - t2.Nodes[i].Start
	}

	depths := make([]int, len(t1.Nodes))
	computeDepth(t1, 0, 1, depths)

	return &CallStructure{
		Quantity: q,
		Sizes1: sizes1, Sizes2: sizes2,
		TwoSidedLeafCalls: make([][]Call, len(t2.Nodes)),
		OneSidedLeafCalls: make([][]Call, len(t2.Nodes)),
		ApproximateCalls: make([][]Call, len(t2.Nodes)),
	}
}

func computeDepth(t *Tree, node, depth int, depths []int) {
	depths[node] = depth
	if t.Nodes[node].Left != -1 {
		computeDepth(t, t.Nodes[node].Left, depth+1, depths)
	}
	if t.Nodes[node].Right != -1 {
		computeDepth(t, t.Nodes[node].Right, depth+1, depths)
	}
}

func (cs *CallStructure) Len() int { return cs.Quantity.Len() }

func (cs *CallStructure) TwoSidedLeaf(t *Tree, i int) {
	cs.TwoSidedLeafCalls[i] = append(cs.TwoSidedLeafCalls[i], Call{ i, -1 })
	cs.Quantity.TwoSidedLeaf(t, i)
}

func (cs *CallStructure) Approximate(t1, t2 *Tree, i2, i1 int) {
	cs.ApproximateCalls[i2] = append(cs.ApproximateCalls[i2], Call{i1, i2})
	cs.Quantity.Approximate(t1, t2, i2, i1)
}

func (cs *CallStructure) OneSidedLeaf(t1, t2 *Tree, i2, i1 int) {
	cs.OneSidedLeafCalls[i2] = append(cs.OneSidedLeafCalls[i2], Call{i1, i2})
	cs.Quantity.OneSidedLeaf(t1, t2, i2, i1)
}
