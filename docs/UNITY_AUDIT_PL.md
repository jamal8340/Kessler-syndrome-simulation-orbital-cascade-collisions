# Audyt dostarczonego projektu Unity

## Zakres audytu

Przeanalizowano zawartość `KesslerSimulation-main.zip`, w
szczególności: - `Assets/SatelliteSpawner.cs`; -
`Assets/OrbitCameraController.cs`; -
`Assets/Scenes/SampleScene.unity`; - `Assets/Satelite.prefab`; -
`Assets/StreamingAssets/satellites_state.json`; -
`Packages/manifest.json`; - `ProjectSettings/ProjectVersion.txt`.

Projekt używa Unity `6000.6.7f1` i URP. Scena zawiera obiekty
`Main Camera`, `Directional Light`, `SatelliteSpawner`, `Global Volume`,
`Manager` i `Earth`.

## Co działa i warto zachować

### Wizualizacja i obsługa

-   Scena 3D z Ziemią i teksturą.
-   Prefab satelity.
-   Kamera orbitalna z obrotem i zoomem.
-   Sterowanie tempem symulacji i pauzą.
-   Wczytywanie tablicy JSON z `StreamingAssets`.
-   Tworzenie obiektów i kolorowanie według kategorii.
-   Panel GUI ze statystykami.
-   Wyłączenie colliderów prefabów, aby Unity PhysX nie próbował
    obsługiwać tysięcy obiektów orbitalnych.

### Ruch orbitalny

`SatelliteSpawner.Update()` przechowuje pozycje i prędkości w km oraz
km/s. Dla każdego obiektu wyznacza przyspieszenie centralnego pola
grawitacyjnego:

`a = -mu * r / |r|^3`

Następnie wykonuje aktualizację prędkości i pozycji:

`v[n+1] = v[n] + a[n] * dt`

`r[n+1] = r[n] + v[n+1] * dt`

To jest metoda symplektycznego Eulera (semi-implicit Euler), a nie
dokładna propagacja SGP4. Jest prostsza i często lepsza energetycznie od
jawnego Eulera, ale pozostaje metodą pierwszego rzędu. Błąd zależy od
kroku czasowego.

## Najważniejsze problemy fizyczne i numeryczne

### 1. Krok czasowy może być za duży

Kod ogranicza delta time klatki, a potem dzieli krok na podkroki
maksymalnie do 10. Przy dużym przyspieszeniu czasu jeden podkrok może
nadal wynosić kilka sekund. Dla szybkich obiektów orbitalnych może to
powodować duży błąd położenia i pominięcie bliskich podejść.

### 2. Hash przestrzenny nie sprawdza sąsiednich komórek

`CheckCollisions()` porównuje obiekt tylko z wcześniejszymi obiektami
przypisanymi do dokładnie tego samego klucza komórki. Nie sprawdza 26
sąsiednich komórek 3D.

Dwa obiekty oddalone o mniej niż próg kolizji mogą leżeć po przeciwnych
stronach granicy komórek i nie zostać porównane. To istotny problem
powodujący fałszywe negatywy.

### 3. Sprawdzana jest odległość chwilowa, nie TCA

Kod sprawdza odległość w chwili wywołania `CheckCollisions()`. Nie
analizuje względnego ruchu między poprzednim a aktualnym krokiem i nie
wyznacza czasu najbliższego podejścia (TCA).

Szybkie spotkanie może nastąpić i zakończyć się pomiędzy kolejnymi
kontrolami.

### 4. Próg 2 km nie jest fizycznym promieniem kolizji

`collisionDistanceKm = 2.0` oznacza, że obiekty w odległości mniejszej
niż 2 km są traktowane jak zderzające się. To może być użyte jako
szeroki próg screeningu koniunkcji, ale nie jako fizyczny warunek
kontaktu. Do fizycznego kontaktu potrzebne są rozmiary obu obiektów i
odpowiednio mała odległość ich środków.

### 5. Filtr minimalnej prędkości względnej

`minRelativeVelocityKmS = 1.0` odrzuca pary o mniejszej prędkości
względnej. Może ograniczać niektóre fałszywe alarmy, ale nie jest
ogólnym warunkiem fizycznym kolizji. Wolniejsze spotkanie również może
zakończyć się kontaktem.

### 6. Obsługa kolizji nie modeluje fragmentacji

`HandleCollision()` oznacza dwa istniejące obiekty jako odłamki.
`DestroySatelliteIntoDebris()` dodaje losowy impuls prędkości do
istniejącego obiektu.

Kod nie: - tworzy nowych fragmentów; - rozdziela masy; - oblicza energii
zderzenia; - zachowuje pędu układu; - wykorzystuje empirycznego modelu
fragmentacji; - propaguje rozkładu rozmiarów i prędkości fragmentów.

Jest to wizualny placeholder zdarzenia, nie fizyczny model kaskady.

### 7. Test klawiszem C jest sztucznym zdarzeniem

`TriggerTestCollision()` wybiera losowy aktywny obiekt w zadanym
przedziale LEO i zmienia jego status. Nie wykrywa rzeczywistego
spotkania dwóch obiektów. Jest przydatny do demonstracji interfejsu, ale
nie może być używany jako wynik naukowy.

### 8. Niezgodność układów odniesienia i osi wymaga sprawdzenia

Pythonowy SGP4 zwraca wektory w układzie TEME. Unity kopiuje składowe
`x, y, z` bez jawnej transformacji. W Unity konwencja osi świata często
używa Y jako osi pionowej, natomiast dane orbitalne zwykle definiują oś
Z jako oś bieguna w przyjętym układzie geocentrycznym.

Należy sprawdzić orientację modelu Ziemi i jawnie zdefiniować
transformację osi. Nie wolno zakładać, że bezpośrednie przypisanie XYZ
jest poprawne dla wizualizacji geograficznej.

Obrót Ziemi wokół `Vector3.up` sam w sobie nie zapewnia poprawnej
transformacji TEME do układu związanego z Ziemią (ECEF).

### 9. Stan wejściowy i epoka

Unity wczytuje statyczny `satellites_state.json`, ale w strukturze
danych używanej przez `SatelliteSpawner` nie ma znacznika epoki ani
informacji o układzie odniesienia. To utrudnia powtarzalność i
porównanie z Pythonem.

### 10. Dane w ZIP-ie różnią się od aktualnego wyniku Pythona

Plik `Assets/StreamingAssets/satellites_state.json` zawiera 18 579
obiektów: - 15 897 `active`; - 2 677 `debris`; - 5 `rocket_body`.

Ostatnie uruchomienie Pythona raportowało 18 479 obiektów. To różne
migawki lub różne etapy przetwarzania. Należy wygenerować jeden
kanoniczny plik i używać tej samej wersji danych w obu systemach.

### 11. Precyzja typów

Unity używa `Vector3` i `float` do przechowywania pozycji w kilometrach.
Dla skali orbity LEO może to być wystarczające do wizualizacji, ale dla
dokładnych obliczeń bliskich podejść lepsze są obliczenia w `double` i
odjęcie wektorów przed konwersją do float.

## Ocena ogólna

  -----------------------------------------------------------------------
  Obszar                              Ocena
  ----------------------------------- -----------------------------------
  Scena 3D i wizualizacja             Dobry fundament prototypu

  Ładowanie katalogu JSON             Działa, ale brakuje metadanych
                                      epoki/układu

  Ruch w centralnym polu              Uproszczony model dwuciałowy
  grawitacyjnym                       

  Integracja numeryczna               Semi-implicit Euler, wymaga testów
                                      kroku

  Screening kolizji                   Niekompletny, ma istotne fałszywe
                                      negatywy

  Fizyczny kontakt                    Niezaimplementowany

  Fragmentacja                        Placeholder, bez fizycznego bilansu

  Kaskada Kesslera                    Demonstracja zdarzeń, nie
                                      zwalidowany model

  Integracja z Pythonem               Jeszcze nie istnieje
  -----------------------------------------------------------------------

## Zalecany kierunek

Nie przepisywać całego Unity. Zachować: - scenę, Ziemię, tekstury i
prefab; - kamerę, sterowanie czasem i GUI; - mechanizm renderowania
obiektów; - wizualne oznaczanie zdarzeń.

Stopniowo zastępować lub podporządkowywać warstwę obliczeniową wynikom
backendu Python: 1. Wspólny format danych z NORAD ID, epoch UTC, frame i
jednostkami. 2. Kanoniczny snapshot obiektów. 3. Pythonowy moduł
close-approach/TCA. 4. Zdarzenia i trajektorie eksportowane do Unity. 5.
Dopiero potem model fragmentacji i kaskady.

Unity może nadal mieć uproszczony tryb demonstracyjny, ale należy go
wyraźnie oznaczyć jako tryb wizualny, nie naukowy.
