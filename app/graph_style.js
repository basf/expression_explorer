[

   {"selector":"node", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "background-color": "white",
       "border-color": "black",
       "content": "data(id)",
       "border-width": "1px",
       "height": "60px",
       "width": "60px"
      }},

   {"selector":"node[lfc<=0]", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "background-color": "mapData(lfc, -3, 0, red, white)",
       "border-color": "black",
       "content": "data(id)",
       "border-width": "1px"
       }},

   {"selector":"node[lfc>0]", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "background-color": "mapData(lfc, 0, 3, white, lightgreen)",
       "border-color": "black",
       "content": "data(id)",
       "border-width": "1px"
       }},

   {"selector":"node:selected", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "border-color": "black",
       "content": "data(id)",
       "border-width": "3px",
       "overlay-opacity": 0.5,
       "overlay-color": "blue"
       }},

    {"selector": "edge", "css": {
        "line-color": "black",
        "source-arrow-shape": "none",
        "target-arrow-shape": "none",
        
        "curve-style": "bezier"
        }},

	{"selector": "edge[type='activation']", "css": {
        "line-color": "gray",
        "source-arrow-shape": "none",
        "target-arrow-shape": "triangle",
        
        "curve-style": "bezier"
        }},
        
     {"selector": "edge[type='repression']", "css": {
        "line-color": "gray",
        "source-arrow-shape": "none",
        "target-arrow-shape": "tee",
        
        "curve-style": "bezier"
        }},

    {"selector": "edge[interaction<0]", "css": {
        "line-color": "mapData(interaction, -1, 0, red, lightGray)"
        }},

    {"selector": "edge[interaction>0]", "css": {
        "line-color": "mapData(interaction, 0, 1, lightGray, green)"
        }},

    {"selector": "edge:selected", "css": {
       "overlay-opacity": 0.2,
       "overlay-color": "maroon"
        }}

   ]
