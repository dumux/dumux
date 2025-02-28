# Custom input data

If you want to read custom user data, let's say to read a couple of check point for your time loop from file, here is some recipes of how to achieve that with simple helpers available in DuMu<sup>x</sup>.

## Read an array of numbers
If you have a simple text file containing a list of numbers
```
0.1
0.2
0.3
0.4
```
You can read these numbers into a `std::vector<double>` like this
```cpp
#include <dumux/io/container.hh>

const auto numbers = Dumux::readFileToContainer<std::vector<double>>("numbers.txt");
```
The file extension is arbitrary. The container type `std::vector<>` needs to support `begin()`, `end()`, and `push_back()`. It has to have an alias `value_type` (the type of the element, here `double`) and this type needs to have an overloaded `operator>>` to be constructible from an input stream. With similar code you can also multi-column data like this
```
0.1 0.2 0.3
0.1 0.4 0.5
```
for example into a `std::vector<Dune::FieldVector<float, 3>>`
```cpp
#include <dune/common/fvector.hh>
#include <dumux/io/container.hh>

const auto coords = Dumux::readFileToContainer<std::vector<Dune::FieldVector<double, 3>>>("coords.txt");
```

## Read INI-style config data
You can read config data (key-value data) using the `Dune::ParameterTreeParser`. This is the same type of file we use for DuMu<sup>x</sup> `input` files. Given a config file like this
```ini
[MyData]
InjectionRate = 0.1 0.2 0.3 0.4
BoundaryFlags = 3 4 5 6

[Problem]
EnableGravity = false
```
the data can be read into a `Dune::ParameterTree` like this
```cpp
#include <dune/common/parametertreeparser.hh>

Dune::ParameterTree config;
Dune::ParameterTreeParser::readINITree("config.ini", config);
```
Again, the file extension is arbitrary. You can retrieve data by casting it into another type like this
```cpp
const bool enableGravity = config.get<bool>("Problem.EnableGravity");
const auto boundaryFlags = config.get<std::vector<int>>("MyData.BoundaryFlags");
```

## Read XML-formatted config data
[XML (Extensible Markup Language)](https://en.wikipedia.org/wiki/XML) is a very flexible structured data format. You can read and write XML files using the [TinyXML-2](http://www.grinninglizard.com/tinyxml2/) library contained in DuMu<sup>x</sup>. For example, given an XML file like this
```xml
<?xml version="1.0"?>
<!-- Simple custom data file -->
<InputData>
  <!-- Injection rates are in kg/s -->
  <InjectionRates>
  0.1 0.2 0.3 0.4
  </InjectionRates>

  <!-- BoundaryTypes: -1: no boundary, 0: Dirichlet, 1: Noflow, 2: Infiltration, 3: Outflow -->
  <BoundaryTypes>
  0 -1 -1 2
  </BoundaryTypes>
</InputData>
```
You can initialise the file reader like this
```cpp
#include <dumux/io/xml/tinyxml2.h>

tinyxml2::XMLDocument xmlData;
const auto returnCode = xmlData.LoadFile("mydata.xml"); // extension arbitrary
if (returnCode != tinyxml2::XML_SUCCESS)
    DUNE_THROW(Dune::IOError, "Couldn't open XML file.");
```
and read fields like this

```cpp
const tinyxml2::XMLElement* inputData = xmlData.FirstChildElement("InputData");
std::stringstream injectionData(inputData->FirstChildElement("InjectionRates")->GetText());
const auto injectionRates = readStreamToContainer<std::vector<double>>(injectionData);
```
