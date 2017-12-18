#make the modules page default view clearer (toggleLevel(1))
sed -i 's/\(init_search();\)/\1 toggleLevel(1);/' html/modules.html
