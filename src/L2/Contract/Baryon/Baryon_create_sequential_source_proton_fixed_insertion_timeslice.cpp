#include "Baryon.ih"

namespace Contract
{
  // this is another version of sequential source, with fixed operator timeslice.
  // It is basically the forward propagator times the operator
  // and the sequential propagator completes the line to the sink which is therefore
  // variable while one inversion is needed for each operator timeslice.
  void create_sequential_source_proton_fixed_insertion_timeslice(Core::Propagator *seqSrc_u,
                                                                 Core::Propagator *seqSrc_d,
                                                                 Core::Propagator const &u,
                                                                 Core::Propagator const &d,
                                                                 size_t const t_op, Base::Operator op)
  {

    size_t const L(u.L());
    size_t const T(u.T());
    assert(L == d.L() && T == d.T());
    assert(L == seqSrc_u->L() && T == seqSrc_u->T());
    assert(L == seqSrc_d->L() && T == seqSrc_d->T());

    Base::Weave weave(L, T);

     (*seqSrc_u) = Core::Propagator(u);
     (*seqSrc_d) = Core::Propagator(d);

     (*seqSrc_u).select_timeslice(t_op);
     (*seqSrc_d).select_timeslice(t_op);

    //if op != 044 

    seqSrc_u->rightMultiplyOperator(op);
    seqSrc_d->rightMultiplyOperator(op);

    
   

  }

}

