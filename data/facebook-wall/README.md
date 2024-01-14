# Facebook wall posts network, part of the Koblenz Network Collection

This directory contains the TSV and related files of the facebook-wosn-wall
network: The is the directed network of a small subset of posts to other user's
wall on Facebook. The nodes of the network are Facebook users, and each directed
edge represents one post, linking the users writing a post to the users whose
wall the post is written on. Since users may write multiple posts on a wall, the
network allows multiple edges connecting a single node pair. Since users may
write on their own wall, the network contains loops.

More information about the network is provided here:
<http://konect.cc/networks/facebook-wosn-wall>.

## Files

- `meta.yaml`: metadata about the network 
- `graph.csv`: the edge list of the network in comma-separated values format,
  with one edge per line

The meaning of the columns in `graph.tsv` is
- `from`: ID of *from* node 
- `to`: ID of *to* node
- `time`: timestamp of edges (UNIX time)

## References

```bib
@MISC{konect:2017:facebook-wosn-wall,
    title = {Facebook wall posts network dataset -- {KONECT}},
    month = oct,
    year = {2017},
    url = {http://konect.cc/networks/facebook-wosn-wall}
}

@inproceedings{viswanath09,
    author = {Viswanath, Bimal and Mislove, Alan and Cha, Meeyoung
              and Gummadi, Krishna P.},
    title = {On the Evolution of User Interaction in {Facebook}},
    booktitle = {Proc. Workshop on Online Soc. Netw.},
    year = {2009},
    pages = {37--42},
}

@inproceedings{viswanath09,
    author = {Viswanath, Bimal and Mislove, Alan and Cha, Meeyoung
              and Gummadi, Krishna P.},
    title = {On the Evolution of User Interaction in {Facebook}},
    booktitle = {Proc. Workshop on Online Soc. Netw.},
    year = {2009},
    pages = {37--42},
}

@inproceedings{konect,
    title = {{KONECT} -- {The} {Koblenz} {Network} {Collection}},
    author = {Jérôme Kunegis},
    year = {2013},
    booktitle = {Proc. Int. Conf. on World Wide Web Companion},
    pages = {1343--1350},
    url = {http://dl.acm.org/citation.cfm?id=2488173},
    url_presentation = {https://www.slideshare.net/kunegis/presentationwow},
    url_web = {http://konect.cc/},
    url_citations = {https://scholar.google.com/scholar?cites=7174338004474749050},
}

@inproceedings{konect,
    title = {{KONECT} -- {The} {Koblenz} {Network} {Collection}},
    author = {Jérôme Kunegis},
    year = {2013},
    booktitle = {Proc. Int. Conf. on World Wide Web Companion},
    pages = {1343--1350},
    url = {http://dl.acm.org/citation.cfm?id=2488173},
    url_presentation = {https://www.slideshare.net/kunegis/presentationwow},
    url_web = {http://konect.cc/},
    url_citations = {https://scholar.google.com/scholar?cites=7174338004474749050},
}
```
