# Kessler Syndrome Simulation -- dokumentacja projektu

> **Status:** prototyp badawczy w rozwoju\
> **Główne technologie:** Python, SGP4, NumPy (planowane), Matplotlib,
> Unity 6\
> **Cel:** analiza populacji obiektów orbitalnych, bliskich podejść,
> kolizji i kaskadowego powstawania odłamków.

## 0. Co nowego w obecnym etapie

W ostatnim etapie rozwoju: - zakończono reorganizację katalogów projektu
Python; - poprawiono ścieżki w pipeline pobierania i przetwarzania
TLE; - dodano walidator TLE z raportami JSON/CSV; - poprawiono parser,
aby pomijał puste linie w plikach CelesTrak; - zaimplementowano
obliczanie elementów oskulacyjnych; - przetworzono 18 479 obiektów bez
pominięć w ostatnim uruchomieniu; - wygenerowano mapę populacji
orbitalnej; - przeprowadzono przegląd istniejącego kodu Unity i
zidentyfikowano różnice między prototypem wizualnym a modelem naukowym.

## 1. Cel projektu

Projekt ma rozwijać symulację zjawiska zespołu Kesslera: wzrostu liczby
obiektów orbitalnych i ryzyka kolejnych zderzeń wskutek generowania
odłamków podczas kolizji.

Docelowo projekt łączy: - rzeczywiste dane orbitalne z katalogów TLE; -
propagację orbit modelem SGP4; - obliczenia mechaniki orbitalnej; -
analizę bliskich podejść i wykrywanie potencjalnych kolizji; - model
fizyki zderzeń i fragmentacji; - długoterminową ewolucję populacji; -
wizualizację 3D w Unity.

**Ważne:** aktualny projekt jest prototypem. Poprawne uruchomienie kodu
i zgodność prostych testów nie oznaczają jeszcze, że cały model jest
zwalidowaną symulacją naukową.

## 2. Architektura: Python i Unity

### Python -- zaplecze obliczeniowe

Python jest rozwijany jako warstwa obliczeniowa odpowiedzialna za: -
pobieranie i przechowywanie danych TLE; - walidację składni i sum
kontrolnych TLE; - propagację stanów orbitalnych SGP4; - wyznaczanie
elementów oskulacyjnych z wektorów stanu; - analizę populacji; -
przyszłe wykrywanie bliskich podejść, obliczanie TCA i modelowanie
kolizji.

### Unity -- istniejąca wizualizacja i prototyp dynamiki

Unity nie jest wyrzucane ani zastępowane. Dostarczony projekt Unity
zawiera już: - scenę z Ziemią, kamerą i obiektem `SatelliteSpawner`; -
teksturę Ziemi i prefab satelity; - wczytywanie `satellites_state.json`
z `StreamingAssets`; - tworzenie obiektów dla katalogu satelitów; -
kolorowanie obiektów według kategorii; - prosty model ruchu w polu
grawitacyjnym Ziemi; - przyspieszanie i pauzowanie czasu; - wstępny
mechanizm wyszukiwania bliskich obiektów; - testowe wywołanie zdarzenia
klawiszem `C`; - prosty interfejs ekranowy ze statystykami.

Obecnie Python i Unity **nie są jeszcze zintegrowane jako jeden wspólny
silnik**. Unity korzysta z własnego pliku JSON i samodzielnie propaguje
ruch obiektów. Docelowo Python ma dostarczać zweryfikowane stany i
zdarzenia, a Unity ma je wizualizować.

## 3. Co zostało wykonane w Pythonie

### 3.1. Struktura i pipeline danych

-   Uporządkowano katalogi na dane surowe, przetworzone, metadane, kod,
    wyniki i dokumentację.
-   Skrypty pobierania i przetwarzania TLE przeniesiono do
    `src/data_pipeline/`.
-   Poprawiono ścieżki tak, aby były wyznaczane względem katalogu
    głównego projektu, a nie położenia skryptu.
-   Dane wejściowe są przechowywane w `data/raw/tle/`.

### 3.2. Pobieranie i walidacja TLE

Zaimplementowano: - pobieranie grup CelesTrak: aktywne obiekty oraz
wybrane grupy odłamków; - zapisywanie grup do oddzielnych plików TLE; -
walidację długości linii i sum kontrolnych; - sprawdzanie zgodności
identyfikatorów NORAD w parach linii; - parsowanie przez bibliotekę
`sgp4`; - test propagacji na epoce TLE; - raport JSON i CSV.

Ostatni przekazany wynik walidacji: \| Metryka \| Wynik \| \|---\|---:\|
\| Pliki TLE \| 5 \| \| Rozpoznane rekordy \| 18 582 \| \| Poprawne
rekordy \| 18 582 \| \| Niepoprawne rekordy \| 0 \| \| Niekompletne
rekordy \| 0 \| \| Ostrzeżenia \| 6 \| \| Duplikaty NORAD między plikami
\| 0 \|

Sześć ostrzeżeń pozostawiono na później. Wynik dotyczy konkretnego
zestawu danych i czasu wykonania walidacji, a nie gwarantuje aktualności
TLE w przyszłości.

### 3.3. Przetwarzanie SGP4

Skrypt `process_TLE_to_json.py` przetwarza rekordy TLE przez SGP4 i
zapisuje pozycję oraz prędkość obiektów w JSON. Są to wektory stanu w
układzie TEME, w kilometrach i kilometrach na sekundę.

Ostatnie raportowane uruchomienie modułu mechaniki orbitalnej wczytało
18 479 obiektów. Rozbieżność względem 18 582 rekordów z raportu
walidacji nie została jeszcze wyjaśniona.

### 3.4. Mechanika orbitalna

W `src/orbital_mechanics/` zaimplementowano operacje wektorowe i
obliczenia: - momentu pędu właściwego; - wektora mimośrodu i
mimośrodu; - energii właściwej; - półosi wielkiej; - inklinacji; -
RAAN; - argumentu perygeum; - anomalii prawdziwej; - wysokości perygeum
i apogeum; - okresu orbitalnego.

Elementy są wyznaczane jako **chwilowe elementy oskulacyjne modelu
dwuciałowego**, obliczone z pozycji i prędkości uzyskanych przez SGP4.
Nie są to niezmienne elementy orbity ani bezpośrednio średnie elementy
TLE.

Test idealnej orbity kołowej: - pozycja: `(7000, 0, 0)` km; - prędkość:
około `(0, 7.54605329, 0)` km/s; - półoś wielka: około 7000 km; -
mimośród: praktycznie 0; - inklinacja: 0°; - perygeum i apogeum: około
629 km; - okres: około 97.14 min.

Test dodatkowej orbity nachylonej i eliptycznej wykazał niezerowy
mimośród i inklinację, ale wybrane warunki początkowe dały perygeum
poniżej powierzchni Ziemi. Był to test matematyczny, a nie fizycznie
dopuszczalna orbita satelity.

Ostatnie raportowane przetwarzanie: \| Metryka \| Wynik \| \|---\|---:\|
\| Obiekty wejściowe \| 18 479 \| \| Przetworzone \| 18 479 \| \|
Pominięte \| 0 \|

Wygenerowane wyniki: - `data/processed/orbital_elements.json` -
`results/figures/orbital_population_map.png`

Wykres pokazuje rozkład wysokości półosi wielkiej względem inklinacji.
Nie jest mapą chwilowych pozycji ani gęstości przestrzennej.

## 4. Audyt dostarczonego projektu Unity

Szczegółowe omówienie znajduje się w
[docs/UNITY_AUDIT_PL.md](docs/UNITY_AUDIT_PL.md).

Najważniejszy wniosek: Unity jest wartościowym prototypem wizualizacji i
uproszczonej dynamiki orbitalnej, ale obecna implementacja kolizji nie
jest jeszcze wiarygodnym modelem fizycznym.

Główne ograniczenia: - ruch jest propagowany niezależnie w Unity metodą
symplektycznego Eulera pierwszego rzędu; - duży krok czasu może
powodować znaczne błędy; - wyszukiwanie kolizji sprawdza tylko obiekty w
tej samej komórce siatki przestrzennej, bez sąsiednich komórek; -
wykrywanie sprawdza chwilową odległość, a nie minimum odległości między
kolejnymi krokami; - próg 2 km oznacza bliskie podejście/koniunkcję, a
nie fizyczny kontakt dwóch obiektów; - obsługa zdarzenia zmienia
etykiety i dodaje losowy impuls prędkości, ale nie generuje nowych
fragmentów i nie zachowuje pędu ani energii; - test klawiszem `C` jest
sztucznym zdarzeniem testowym, a nie wynikiem obliczenia spotkania
orbitalnego; - plik JSON dołączony do Unity zawiera 18 579 obiektów, co
różni się od bieżącego wyniku Pythona 18 479; dane należy
zsynchronizować; - obecnie nie ma integracji Python--Unity.

Nie oznacza to, że Unity jest bezużyteczne. Scena, prefabrykaty, kamera,
prezentacja populacji, interfejs i mechanizm animacji są dobrym
fundamentem.

## 5. Ograniczenia naukowe obecnego modelu

1.  SGP4 zwraca stany w układzie TEME. Przed porównywaniem z innymi
    danymi lub użyciem w modelu Ziemi należy jawnie ustalić układ
    odniesienia, epokę i transformacje.
2.  Pozycje i prędkości w pliku JSON muszą mieć wspólny czas
    odniesienia.
3.  Elementy oskulacyjne wyliczane z chwilowego stanu nie powinny być
    utożsamiane ze średnimi elementami TLE.
4.  Model dwuciałowy nie uwzględnia automatycznie J2, oporu
    atmosferycznego, ciśnienia promieniowania ani perturbacji od
    Księżyca i Słońca.
5.  Obecne masy i promienie przypisywane kategoriom obiektów są
    wartościami szacunkowymi, nie pomiarami poszczególnych satelitów.
6.  Nie ma jeszcze obliczeń prawdopodobieństwa kolizji, modelu
    niepewności ani zwalidowanej fragmentacji.

## 6. Plan dalszych prac

### Etap 4 -- spójność danych i analiza populacji

-   Uzgodnić liczbę rekordów w walidacji, JSON Pythona i JSON Unity.
-   Zapisywać czas/epokę danych i układ odniesienia w metadanych.
-   Sprawdzić wykres populacji i rozkłady wysokości oraz inklinacji.
-   Dodać testy jednostkowe dla obliczeń elementów orbitalnych.

### Etap 5 -- Collision Detection

-   Propagować pary obiektów do wspólnego czasu i układu odniesienia.
-   Wprowadzić wydajne filtrowanie kandydatów.
-   Obliczać względną pozycję i prędkość.
-   Wyznaczać czas najbliższego podejścia (TCA) i miss distance.
-   Uwzględnić ciągłość ruchu między krokami, aby nie pomijać szybkich
    spotkań.
-   Rozróżniać bliskie podejście od fizycznego zderzenia.

### Etap 6 -- prawdopodobieństwo kolizji

-   Modelować niepewność pozycji i prędkości.
-   Wprowadzić rozmiary obiektów i geometrię spotkania.
-   Obliczać prawdopodobieństwo kolizji dopiero po zdefiniowaniu modelu
    niepewności.

### Etap 7 -- fizyka zderzeń i fragmentacja

-   Obliczać energię i prędkość względną zderzenia.
-   Wybrać udokumentowany model fragmentacji.
-   Generować rozkład mas, rozmiarów i prędkości odłamków.
-   Kontrolować bilanse masy i pędu w granicach założeń modelu.

### Etap 8 -- propagacja odłamków

-   Propagować nowe fragmenty.
-   Uwzględniać odpowiednie perturbacje i opór atmosferyczny.
-   Aktualizować populację w czasie.

### Etap 9 -- kaskada Kesslera

-   Uruchamiać kolejne zdarzenia kolizyjne.
-   Rejestrować zdarzenia i pochodzenie fragmentów.
-   Porównywać scenariusze i analizować wrażliwość wyników.

### Etap 10 -- integracja z Unity

-   Ustalić jeden format wymiany danych.
-   Dodać stabilne identyfikatory NORAD i znaczniki czasu.
-   Unity powinno wizualizować stany i zdarzenia obliczone przez
    backend, zamiast prowadzić niezależną, rozbieżną symulację.
-   Zachować Unity jako interaktywną warstwę 3D.

## 7. Najbliższy konkretny cel

Najpierw należy: 1. Uzgodnić i zsynchronizować dane wejściowe
Python/Unity. 2. Poprawić i przetestować wykrywanie bliskich podejść w
Pythonie. 3. Zweryfikować wyniki na kontrolowanych przypadkach. 4.
Dopiero potem zasilać Unity rzeczywistymi zdarzeniami koniunkcji i
kolizji.

**Kryterium postępu:** każdy moduł ma mieć opisane założenia, testy i
znane ograniczenia. Sama animacja lub brak błędów wykonania nie stanowią
walidacji naukowej.
