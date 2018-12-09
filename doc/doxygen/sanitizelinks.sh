#make the modules page default view clearer (toggleLevel(1))
if [ -e html/modules.html ]; then
  sed -i 's/\(init_search();\)/\1 toggleLevel(1);/' html/modules.html
fi
if [ -e html/modules.HTML ]; then
  sed -i 's/\(init_search();\)/\1 toggleLevel(1);/' html/modules.HTML
fi
