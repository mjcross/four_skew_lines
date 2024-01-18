from vector import Point, Vector
from fractions import Fraction


def main():

    a0 = Vector(19, 13, 30)
    da = Vector(-2, 1, -2)

    b0 = Vector(18, 19, 22)
    db = Vector(-1, -1, -2)

    c0 = Vector(20, 19, 15)
    dc = Vector(1, -5, -3)

    d0 = Vector(12, 31, 28)
    dd = Vector(-1, -2, -1)

    e0 = Vector(20, 25, 34)
    de = Vector(-2, -2, -4)


    print(f'n = {da.cross(a0 - b0)} - beta * {da.cross(db)}')
    print(f'gamma = ({(a0-c0).dot(da.cross(a0-b0))} - beta * {(a0-c0).dot(da.cross(db))}) / ({dc.dot(da.cross(a0-b0))} - beta * {dc.dot(da.cross(db))})')
    print(f'delta = ({(a0-d0).dot(da.cross(a0-b0))} - beta * {(a0-d0).dot(da.cross(db))}) / ({dd.dot(da.cross(a0-b0))} - beta * {dd.dot(da.cross(db))})')


    for beta in range(10):
        b = b0 + beta * db

        # normal to the plane containing line A and point B
        n = da.cross(a0 - b0) - beta * (da.cross(db))
        assert n.dot(da) == 0
        assert n.dot(a0 - (b0 + beta * db)) == 0

        # intersection of the plane and line C
        if dc.dot(n) == 0:
            continue    # plane is parallel to line C
        else:
            gamma = Fraction((a0 - c0).dot(n), dc.dot(n))
            assert gamma == Fraction((a0-c0).dot(da.cross(a0-b0)) - beta * (a0-c0).dot(da.cross(db)), dc.dot(da.cross(a0-b0)) - beta * dc.dot(da.cross(db)))

            c = c0 + gamma * dc
            assert c.dot(n) == a0.dot(n)

        # intersection of the plane and line 'd'
        if dd.dot(n) == 0:
            continue    # plane is parallel to line D
        else:
            delta = Fraction((a0 - d0).dot(n), dd.dot(n))
            assert delta == Fraction((a0-d0).dot(da.cross(a0-b0)) - beta * (a0-d0).dot(da.cross(db)), dd.dot(da.cross(a0-b0)) - beta * dd.dot(da.cross(db)))
            d = d0 + delta * dd
            assert d.dot(n) == a0.dot(n)

        # check for co-linearity
        bc = c - b
        bd = d - b
        print(beta, bc.isparallel(bd), b.cross(c-d) == d.cross(c)) 



if __name__ == '__main__':
    main()