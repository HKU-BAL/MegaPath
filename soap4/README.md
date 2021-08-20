SOAP4-CPU
===
This is CPU version of SOAP4, which implements exactly the same algorithm in SOAP4.

Main Changes
===
- GPU round 1/2 is removed, as it is unnecessary
- seeds aligned in 0MM will not continue to 1MM anymore
- add varsize seeding 
- add NBM seeding
- add a missing Case in seeding (bidirectional backward lookup)
- add deduplication of NBM only seeds which appear in both Case A and B
- Lookup Table for exact match step in 1MM is removed, as it cannot handle NBM
- fix NBM and mismatch logic, mismatch can be used as an NBM
- fix NBM and mismatch logic, use NBM first when both NBM and mismatch is available
- modify AlignmentResult structure to support varsize seeding and NBM only deduplication
- fix a bug in DP which results in wrong mapping quality
- optimize some SIMD instructions



Known Issues
===
- Check and Extend is removed


### What is this repository for? ###

### How do I get set up? ###

### Contribution guidelines ###

### Who do I talk to? ###