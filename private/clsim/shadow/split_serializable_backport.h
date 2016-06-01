#ifndef SPLIT_SERIALIZABLE_BACKPORT_H_INCLUDED
#define SPLIT_SERIALIZABLE_BACKPORT_H_INCLUDED
// TODO: This file needs to go away once a more recent icetray version makes it
// into a release. It's a backport of "I3_SPLIT_SERIALIZABLE".

#ifndef I3_SPLIT_SERIALIZABLE
#define I3_SPLIT_SERIALIZABLE(T)                                        \
  I3_SERIALIZABLE(T)                                                    \
  template void T::save(icecube::archive::portable_binary_oarchive&, unsigned) const; \
  template void T::load(icecube::archive::portable_binary_iarchive&, unsigned); \
  template void T::load(icecube::archive::xml_iarchive&, unsigned);       \
  template void T::save(icecube::archive::xml_oarchive&, unsigned) const;
#endif

#endif  // SPLIT_SERIALIZABLE_BACKPORT_H_INCLUDED
