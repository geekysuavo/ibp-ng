
/* include the molecular topology header. */
#include "topol.h"

/* include the peptide headers. */
#include "peptide-atoms.h"
#include "peptide-bonds.h"
#include "peptide-angles.h"
#include "peptide-torsions.h"
#include "peptide-impropers.h"

/* topol_apply(): apply a named residue topology entry to a peptide at
 * a specified location in the peptide sequence.
 *
 * arguments:
 *  @top: pointer to the topology structure to access.
 *  @resname: string name of the residue topology to use.
 *  @P: pointer to the peptide structure to modify.
 *  @ires: peptide sequence index.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_apply (topol_t *top, const char *resname,
                 peptide_t *P, unsigned int ires) {
  /* declare required variables:
   *  @res: residue topology structure pointer.
   *  @dihed: dihedral topology structure pointer.
   *  @angle: angle topology structure pointer.
   *  @atom: atom topology structure pointer.
   *  @bond: bond topology structure pointer.
   *  @i: general-purpose loop counter.
   */
  topol_residue_t *res;
  topol_dihedral_t *dihed;
  topol_angle_t *angle;
  topol_atom_t *atom;
  topol_bond_t *bond;
  unsigned int i;

  /* search for the named residue topology entry. */
  for (i = 0; i < top->n_res; i++) {
    /* if a match is found, terminate the search. */
    if (strcmp(top->res[i].name, resname) ==  0)
      break;
  }

  /* check that a valid residue topology entry was found. */
  if (i >= top->n_res)
    throw("unable to locate residue topology for '%s'", resname);

  /* obtain a pointer to the residue topology structure. */
  res = top->res + i;

  /* loop over the atom topologies. */
  for (i = 0; i < res->n_atoms; i++) {
    /* obtain a pointer to the atom topology structure. */
    atom = res->atoms + i;

    /* determine the type of action to perform. */
    switch (atom->mode) {
      /* add a new atom. */
      case TOPOL_MODE_ADD:
        if (!peptide_atom_add(P, ires + atom->off,
                              atom->name, atom->type,
                              atom->mass, atom->charge, 0.0))
          throw("unable to add atom %u.%s (%s)",
                ires + atom->off,
                atom->name, atom->type);

        /* break the switch. */
        break;

      /* modify an existing atom. */
      case TOPOL_MODE_MODIFY:
        if (!peptide_atom_modify(P, ires + atom->off,
                                 atom->name, atom->type,
                                 atom->mass, atom->charge, 0.0))
          throw("unable to modify atom %u.%s (%s)",
                ires + atom->off,
                atom->name, atom->type);

        /* break the switch. */
        break;

      /* delete an existing atom. */
      case TOPOL_MODE_DELETE:
        /* delete all bonds to the atom. */
        if (!peptide_bond_delete_any(P, ires + atom->off, atom->name))
          throw("unable to delete bonds for %u.%s",
                ires + atom->off, atom->name);

        /* delete all angles containing the atom. */
        if (!peptide_angle_delete_any(P, ires + atom->off, atom->name))
          throw("unable to delete angles for %u.%s",
                ires + atom->off, atom->name);

        /* delete all torsions containing the atom. */
        if (!peptide_torsion_delete_any(P, ires + atom->off, atom->name))
          throw("unable to delete dihedrals for %u.%s",
                ires + atom->off, atom->name);

        /* delete all impropers containing the atom. */
        if (!peptide_improper_delete_any(P, ires + atom->off, atom->name))
          throw("unable to delete impropers for %u.%s",
                ires + atom->off, atom->name);

        /* delete the atom. */
        if (!peptide_atom_delete(P, ires + atom->off, atom->name))
          throw("unable to delete atom %u.%s",
                ires + atom->off, atom->name);

        /* break the switch. */
        break;

      /* throw an exception. */
      default:
        throw("unsupported atom topology mode");
    }
  }

  /* loop over the bond topologies. */
  for (i = 0; i < res->n_bonds; i++) {
    /* obtain a pointer to the bond topology structure. */
    bond = res->bonds + i;

    /* determine the type of action to perform. */
    switch (bond->mode) {
      /* add a new bond. */
      case TOPOL_MODE_ADD:
        if (!peptide_bond_add(P, ires + bond->off[0], bond->atoms[0],
                                 ires + bond->off[1], bond->atoms[1], 0))
          throw("unable to add bond (%u.%s -- %u.%s)",
                ires + bond->off[0], bond->atoms[0],
                ires + bond->off[1], bond->atoms[1]);

        /* break the switch. */
        break;

      /* delete an existing bond. */
      case TOPOL_MODE_DELETE:
        if (!peptide_bond_delete(P, ires + bond->off[0], bond->atoms[0],
                                    ires + bond->off[1], bond->atoms[1]))
          throw("unable to delete bond (%u.%s -- %u.%s)",
                ires + bond->off[0], bond->atoms[0],
                ires + bond->off[1], bond->atoms[1]);

        /* break the switch. */
        break;

      /* throw an exception. */
      default:
        throw("unsupported bond topology mode");
    }
  }

  /* loop over the angle topologies. */
  for (i = 0; i < res->n_angles; i++) {
    /* obtain a pointer to the angle topology structure. */
    angle = res->angles + i;

    /* determine the type of action to perform. */
    switch (angle->mode) {
      /* add a new angle. */
      case TOPOL_MODE_ADD:
        if (!peptide_angle_add(P,
              ires + angle->off[0], angle->atoms[0],
              ires + angle->off[1], angle->atoms[1],
              ires + angle->off[2], angle->atoms[2]))
          throw("unable to add angle (%u.%s, %u.%s, %u.%s)",
                ires + angle->off[0], angle->atoms[0],
                ires + angle->off[1], angle->atoms[1],
                ires + angle->off[2], angle->atoms[2]);

        /* break the switch. */
        break;

      /* delete an existing angle. */
      case TOPOL_MODE_DELETE:
        if (!peptide_angle_delete(P,
              ires + angle->off[0], angle->atoms[0],
              ires + angle->off[1], angle->atoms[1],
              ires + angle->off[2], angle->atoms[2]))
          throw("unable to delete angle (%u.%s, %u.%s, %u.%s)",
                ires + angle->off[0], angle->atoms[0],
                ires + angle->off[1], angle->atoms[1],
                ires + angle->off[2], angle->atoms[2]);

        /* break the switch. */
        break;

      /* throw an exception. */
      default:
        throw("unsupported angle topology mode");
    }

  }

  /* loop over the torsion topologies. */
  for (i = 0; i < res->n_torsions; i++) {
    /* obtain a pointer to the torsion topology structure. */
    dihed = res->torsions + i;

    /* determine the type of action to perform. */
    switch (dihed->mode) {
      /* add a new torsion. */
      case TOPOL_MODE_ADD:
        if (!peptide_torsion_add(P,
              ires + dihed->off[0], dihed->atoms[0],
              ires + dihed->off[1], dihed->atoms[1],
              ires + dihed->off[2], dihed->atoms[2],
              ires + dihed->off[3], dihed->atoms[3]))
          throw("unable to add dihedral (%u.%s, %u.%s, %u.%s, %u.%s)",
                ires + dihed->off[0], dihed->atoms[0],
                ires + dihed->off[1], dihed->atoms[1],
                ires + dihed->off[2], dihed->atoms[2],
                ires + dihed->off[3], dihed->atoms[3]);

        /* break the switch. */
        break;

      /* delete an existing torsion. */
      case TOPOL_MODE_DELETE:
        if (!peptide_torsion_delete(P,
              ires + dihed->off[0], dihed->atoms[0],
              ires + dihed->off[1], dihed->atoms[1],
              ires + dihed->off[2], dihed->atoms[2],
              ires + dihed->off[3], dihed->atoms[3]))
          throw("unable to add dihedral (%u.%s, %u.%s, %u.%s, %u.%s)",
                ires + dihed->off[0], dihed->atoms[0],
                ires + dihed->off[1], dihed->atoms[1],
                ires + dihed->off[2], dihed->atoms[2],
                ires + dihed->off[3], dihed->atoms[3]);

        /* break the switch. */
        break;

      /* throw an exception. */
      default:
        throw("unsupported dihedral topology mode");
    }
  }

  /* loop over the improper topologies. */
  for (i = 0; i < res->n_impropers; i++) {
    /* obtain a pointer to the improper topology structure. */
    dihed = res->impropers + i;

    /* determine the type of action to perform. */
    switch (dihed->mode) {
      /* add a new improper. */
      case TOPOL_MODE_ADD:
        if (!peptide_improper_add(P,
              ires + dihed->off[0], dihed->atoms[0],
              ires + dihed->off[1], dihed->atoms[1],
              ires + dihed->off[2], dihed->atoms[2],
              ires + dihed->off[3], dihed->atoms[3]))
          throw("unable to add improper (%u.%s, %u.%s, %u.%s, %u.%s)",
                ires + dihed->off[0], dihed->atoms[0],
                ires + dihed->off[1], dihed->atoms[1],
                ires + dihed->off[2], dihed->atoms[2],
                ires + dihed->off[3], dihed->atoms[3]);

        /* break the switch. */
        break;

      /* delete an existing improper. */
      case TOPOL_MODE_DELETE:
        if (!peptide_improper_delete(P,
              ires + dihed->off[0], dihed->atoms[0],
              ires + dihed->off[1], dihed->atoms[1],
              ires + dihed->off[2], dihed->atoms[2],
              ires + dihed->off[3], dihed->atoms[3]))
          throw("unable to add improper (%u.%s, %u.%s, %u.%s, %u.%s)",
                ires + dihed->off[0], dihed->atoms[0],
                ires + dihed->off[1], dihed->atoms[1],
                ires + dihed->off[2], dihed->atoms[2],
                ires + dihed->off[3], dihed->atoms[3]);

        /* break the switch. */
        break;

      /* throw an exception. */
      default:
        throw("unsupported improper topology mode");
    }
  }

  /* return success. */
  return 1;
}

/* topol_apply_all(): add topology information to a peptide structure.
 *
 * arguments:
 *  @top: pointer to the topology structure to access.
 *  @P: pointer to the peptide structure to modify.
 *
 * returns:
 *  integer indicating whether (1) or not (0) the operation succeeded.
 */
int topol_apply_all (topol_t *top, peptide_t *P) {
  /* declare required variables:
   *  @i: general-purpose loop counter.
   *  @ret: return status variable.
   */
  unsigned int i;
  int ret;

  /* loop over the peptide sequence to make some sidechains explicit. */
  for (i = 0; i < P->n_res; i++) {
    /* if the current residue is proline, make it explicit. */
    const char *resi = peptide_get_resname(P, i);
    if (strcmp(resi, "PRO") == 0 && !peptide_add_sidechain(P, i))
      throw("failed to make residue %u (%s) explicit", i + 1, resi);
  }

  /* loop over the peptide sequence to add each residue topology. */
  for (i = 0; i < P->n_res; i++) {
    /* add the currently indexed residue topology. */
    if (!topol_apply(top, peptide_get_restype(P, i), P, i))
      throw("unable to apply topology to %s%u",
            peptide_get_resname(P, i), i + 1);
  }

  /* loop over the peptide sequence to link adjacent residues. */
  for (i = 0; i < P->n_res - 1; i++) {
    /* link the current residue to its next neighbor. */
    if (strcmp(peptide_get_resname(P, i + 1), "PRO") == 0) {
      /* use a linkage to proline. */
      ret = topol_apply(top, "PEPP", P, i);
    }
    else {
      /* use a standard linkage. */
      ret = topol_apply(top, "PEPT", P, i);
    }

    /* check for errors. */
    if (!ret)
      throw("unable to link %s%u to %s%u",
            peptide_get_resname(P, i), i + 1,
            peptide_get_resname(P, i + 1), i + 2);
  }

  /* apply n-terminal patches to the first residue. */
  if (!topol_apply(top, "NTER", P, 0))
    throw("unable to patch n-terminal topology");

  /* apply c-terminal patches to the last residue. */
  if (!topol_apply(top, "CTER", P, P->n_res - 1))
    throw("unable to patch c-terminal topology");

  /* return success. */
  return 1;
}

