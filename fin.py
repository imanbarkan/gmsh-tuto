import argparse
import gmsh
import sys


def main():

    MeshAlgo2D = {
        "MeshAdapt": 1,
        "Automatic": 2,
        "Initial": 3,
        "Delaunay": 5,
        "Frontal-Delaunay": 6,
        "BAMG": 7,
    }

    """Console script for python_magnetgeo."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lc", help="load mesh size from file", type=float, default=0.2)
    parser.add_argument("--mesh", help="activate mesh mode",
                        action="store_true")
    parser.add_argument(
        "--fragment", help="activate fragment", action="store_true")
    parser.add_argument("--view", help="activate view mode",
                        action="store_true")
    parser.add_argument(
        "--debug", help="activate debug mode", action="store_true")

    args = parser.parse_args()
    lc = args.lc

    # init gmsh
    gmsh.initialize()

    # add a model
    gmsh.model.add("fin")

    Nfins = 4
    Lfins = 2.5
    t = 0.25
    d = 0.75

    x0 = 0
    y0 = 0
    z0 = 0
    dx = 1
    dy = 0.25
    _id = gmsh.model.occ.addRectangle(x0, y0, z0, dx, Nfins * (d+t) + t)
    gmsh.model.occ.synchronize()

    # create a branch

    id_list = []
    for i in range(1, Nfins + 1):
        _id_branch = gmsh.model.occ.addRectangle(
            x0-Lfins, y0 + i * (d+t), z0, 2*Lfins + dx, t)
        id_list.append(_id_branch)
        print(f'id_list={id_list}')
    # _id_branch = gmsh.model.occ.addRectangle(x0+dx, y0+dy/2., z0, 6*dx, dy/8.)
    # _id_branch = gmsh.model.occ.addRectangle(-x0+dx-dy, y0+dy/2., z0, 6*dx, dy/8.)

    gmsh.model.occ.synchronize()

    if args.fragment:
        ov, ovv = gmsh.model.occ.fragment(
            [(2, _id)], [(2, i) for i in id_list]
        )

        # ov contains all the generated entities of the same dimension as the input
        # entities:
        print("fragment produced volumes:")
        for e in ov:
            print(e)

        # ovv contains the parent-child relationships for all the input entities:
        print("before/after fragment relations:")
        # for e in zip([(2, _id)] + [(2, _id_branch)], ovv):
        for e in zip([(2, _id)] + [(2, i) for i in id_list], ovv):
            print("parent " + str(e[0]) + " -> child " + str(e[1]))

        gmsh.model.occ.synchronize()

    # select 0 for nodes, 1 for lines
    select = 1
    eps = 1.e-3
    xmin = x0 - eps
    ymin = y0 - eps
    zmin = z0 - eps
    xmax = x0 + dx + eps
    ymax = y0 + eps
    zmax = z0 + eps
    _ov = gmsh.model.getEntitiesInBoundingBox(
        xmin, ymin, zmin, xmax, ymax, zmax, select)
    print(f'_ov={_ov}')

    """
    # create a physical group for surface _id
    ps = gmsh.model.addPhysicalGroup(2, [_id])
    gmsh.model.setPhysicalName(2, ps, name="My surface")
    ps = gmsh.model.addPhysicalGroup(2, id_list)
    gmsh.model.setPhysicalName(2, ps, name="My branch")
    """

    surf_princ = []
    for i in range(1, 2*Nfins + 1):
        surf_princ.append(i)
    ps = gmsh.model.addPhysicalGroup(2, surf_princ)
    gmsh.model.setPhysicalName(2, ps, name="Post")

    fin_id = 1
    for i in range(4*Nfins, 2*Nfins, -2):
        ps = gmsh.model.addPhysicalGroup(2, [i, i-1])
        gmsh.model.setPhysicalName(2, ps, name="Fin_"+str(fin_id))
        fin_id += 1

    root = _ov[0][1]

    ps = gmsh.model.addPhysicalGroup(1, [root])
    gmsh.model.setPhysicalName(1, ps, name="Gamma_root")

    ext = gmsh.model.getBoundary(gmsh.model.getEntities(2))
    for i in range(len(ext)):
        ext[i] = ext[i][1]
    ext = [indice for indice in ext if indice != root]

    ps = gmsh.model.addPhysicalGroup(1, ext)
    gmsh.model.setPhysicalName(1, ps, name="Gamma_ext")

    if args.mesh:

        # Frontal-Delaunay for 2D meshes
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.setAlgorithm(2, 1, 1)

        # set mesh characteristic size
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)

        # Affinage de la maille
        gmsh.model.mesh.field.add("Distance", 1)
        # gmsh.model.mesh.field.setNumbers(1, "PointsList", [])
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", ext)
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", lc/10)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lc)
        gmsh.model.mesh.field.setNumber(2, "DistMin", t/4)
        gmsh.model.mesh.field.setNumber(2, "DistMax", t/2)

        gmsh.model.mesh.field.setAsBackgroundMesh(2)

        # set mesh characteristic size
        # gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)
        # We can then generate a 2D mesh...
        gmsh.model.mesh.generate(2)

        # ... and save it to disk
        gmsh.write("fin.msh")

    # to view the model
    if args.view:
        gmsh.fltk.run()

    # need to end gmsh
    gmsh.finalize()

    return 0


if __name__ == "__main__":
    sys.exit(main())
