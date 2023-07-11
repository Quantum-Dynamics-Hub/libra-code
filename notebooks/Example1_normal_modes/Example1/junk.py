si8="""8

Si 0.000000000000000 0.000000000000000 0.000000000000000
Si 0.000000000000000 5.131267690524859 5.131267690524859
Si 5.131267690524859 5.131267690524859 0.000000000000000
Si 5.131267690524859 0.000000000000000 5.131267690524859
Si 7.696901534842424 2.565633844317566 7.696901534842424
Si 2.565633844317566 2.565633844317566 2.565633844317566
Si 2.565633844317566 7.696901534842424 7.696901534842424
Si 7.696901534842424 7.696901534842424 2.565633844317566"""

view1 = py3Dmol.view(width=800,height=400)
view1.addModel(si8,'xyz',{'vibrate': {'frames':10,'amplitude':500}})
view1.setBackgroundColor('0xeeeeee')
view1.setStyle({'sphere':{}})
view1.animate({'loop': 'backAndForth'})
view1.zoomTo()
view1.show()

print si8
print xyz00v


"""
C = QE_methods.read_md_data_cell("x0.xml")

print C.get(0,0), C.get(1,0), C.get(2,0)
print C.get(3,0), C.get(4,0), C.get(5,0)
print C.get(6,0), C.get(7,0), C.get(8,0)

H = MATRIX3x3()
H.xx = C.get(0,0);  H.xy = C.get(0,1);  H.xz = C.get(0,2);
H.yx = C.get(1,0);  H.yy = C.get(1,1);  H.yz = C.get(1,2);
H.zx = C.get(2,0);  H.zy = C.get(2,1);  H.zz = C.get(2,2);

#alat = 10.26253537973
tH = H.T()  #* alat # * units.Angst

fR = fold_coords(R, tH, "abc")
"""