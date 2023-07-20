def matrix_chain_4(m,n,k,l,q):

    A = Matrix("A", (n, m))
    A.set_property(Property.FULL_RANK)

    B = Matrix("B", (m, k))
    B.set_property(Property.FULL_RANK)

    C = Matrix("C", (k, l))
    C.set_property(Property.FULL_RANK)
    
    D = Matrix("D", (l, q))
    D.set_property(Property.FULL_RANK)

    Y = Matrix("Y", (n, q))

    # Y = ABCD
    equations = Equations(Equal(Y, Times(A,B,C,D)))

    return equations