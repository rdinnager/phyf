# format() formats the pf, pfc, and pfp objects correctly

    Code
      format(pfc_test)
    Output
      [1] "(o)--  0.3--> Node2 --  0.1--> A "                 
      [2] "(o)--  0.3--> Node2 --  0.2--> B "                 
      [3] "(o)--  0.8--> Node3 --  0.4--> C "                 
      [4] "(o)--  0.8--> Node3 --  0.7--> Node4 --  0.5--> D "
      [5] "(o)--  0.8--> Node3 --  0.7--> Node4 --  0.6--> E "
      [6] "(o)--  0.3--> Node2 "                              
      [7] "(o)--  0.8--> Node3 "                              
      [8] "(o)--  0.8--> Node3 --  0.7--> Node4 "             

---

    Code
      format(pfp_test)
    Output
      [1] "6->7->1"    "6->7->2"    "6->8->3"    "6->8->9->4" "6->8->9->5"

