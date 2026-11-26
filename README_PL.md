# Kessler Syndrome Simulation – backend

Symulacja zderzeń obiektów na orbicie okołoziemskiej: jak dwa obiekty rozpadają się na fragmenty, co z nich zostaje, dokąd lecą odłamki i jak chmura odłamków rozprzestrzenia się w czasie. Wykorzystuje prawdziwy katalog ok. 18 500 śledzonych satelitów i odłamków (dane CelesTrak) oraz standardowy model rozpadu NASA.

To repozytorium to **backend**: cała fizyka działa tutaj. Wizualizacja 3D to projekt Unity `KesslerSimulation` (patrz jego `README_PL.md`), który łączy się z serwerem tego backendu. Bez Unity z backendu można korzystać z wiersza poleceń.

## Możliwości

- Pokazuje, gdzie znajduje się każdy obiekt z katalogu w dowolnej chwili, w przód i wstecz w czasie.
- Znajduje pary obiektów, które miną się blisko siebie.
- Zderza dowolne dwa obiekty; powstałe fragmenty stają się nowymi obiektami na własnych orbitach (dwie chmury odłamków, po jednej na rodzica); fragmenty, które od razu wpadłyby w atmosferę, są usuwane.
- Cofanie i odtwarzanie: przed zderzeniem pierwotne obiekty wracają; ponowne odtworzenie daje te same fragmenty.
- Odtwarza dwa prawdziwe zdarzenia: zderzenie Iridium 33 × Cosmos 2251 (2009) i test broni antysatelitarnej na Fengyun‑1C (2007), oraz porównuje symulowane odłamki z prawdziwymi.
- Udostępnia to wszystko na żywo wizualizacji w Unity.

## Wymagania i instalacja

- [uv](https://docs.astral.sh/uv/) – przy pierwszym `uv run` samo instaluje właściwą wersję Pythona (3.14) i wszystkie zależności.
- Dane potrzebne do działania są w repozytorium (`data/`).

```bash
uv sync          # opcjonalnie: instaluje wszystko z góry
```

Wszystkie poniższe polecenia uruchamia się w folderze tego repozytorium.

## Długoterminowy opór atmosferyczny (F1, eksperymentalny)

Opcjonalny model oporu uśrednionego po orbicie zachowuje historię orbit i momenty re-entry. Nie jest włączony w istniejącym serwerze live; przeznaczono go dla przyszłego trybu wsadowego kaskady. Uruchom dwa fragmenty referencyjne lub pełny katalog bez kolizji:

```bash
uv run python -m kessler_sim.drag_demo --years 50
uv run python -m kessler_sim.drag_demo --scenario catalogue --years 10
uv run python -m kessler_sim.drag_demo --scenario catalogue --years 10 --no-stationkeeping
uv run pytest tests/test_drag.py tests/test_stationkeeping.py -q
```

Domyślnie nieuszkodzone aktywne satelity idealnie kompensują opór (stała półoś wielka i mimośród); odłamki, śmieci i pozostałości po kolizjach pozostają pasywne. Użyj `--no-stationkeeping` do porównania z pierwotnym modelem pasywnym. To założenie bezterminowego działania, a nie model rzeczywistych czasów misji.

Założenia, API Pythona i ograniczenia opisano w `docs/ATMOSPHERIC_DRAG.md`.

## Szybki start z Unity

```bash
uv run python -m kessler_sim.server --play --speed 600
```

Następnie otwórz projekt Unity i naciśnij Play (patrz `README_PL.md` w Unity). Gdy serwer jest gotowy, wypisuje `listening on ws://127.0.0.1:8765`; prowadzi jedną symulację wspólną dla wszystkich podłączonych. Zatrzymuje go Ctrl+C; aby zmienić scenariusz albo zacząć od nowa, uruchom go ponownie.

### Opcje serwera

| Opcja | Domyślnie | Znaczenie |
|---|---|---|
| `--scenario NAZWA` | `catalogue` | Który scenariusz wczytać (patrz niżej) |
| `--play` | wyłączone | Start z uruchomionym zegarem |
| `--speed S` | 60 | Sekundy symulacji na sekundę rzeczywistą |
| `--port P` | 8765 | Port, z którym łączy się Unity |
| `--host H` | 127.0.0.1 | Adres sieciowy, na którym serwer nasłuchuje |
| `--rate R` | 20 | Liczba aktualizacji pozycji na sekundę wysyłanych do Unity |

## Scenariusze

| Scenariusz | Co pokazuje |
|---|---|
| `catalogue` | Wszystkie 18 558 obiektów katalogu w chwili z zestawu danych (2026‑11‑12 22:50 UTC); bez zderzeń – pole do własnych zderzeń |
| `closest_approach_breakup` | Pełny katalog; najbliższa para przelatująca w ciągu pierwszych 10 minut zderza się (obecnie STARLINK‑6135 × SIGI, ok. 766 fragmentów) |
| `iridium_cosmos_2009` | Odtworzenie zderzenia Iridium 33 i Cosmos 2251 z 2009 r. (tylko te dwa satelity i ich odłamki); zderzenie w t = 600 s |
| `fengyun_1c_2007` | Odtworzenie testu broni antysatelitarnej z 2007 r.: pocisk przechwytujący wznosi się z Ziemi i niszczy satelitę meteorologicznego Fengyun‑1C; zderzenie w t = 600 s |

`uv run python -m kessler_sim.demo scenarios` wypisuje ich listę. Scenariusze historyczne używają wyidealizowanych orbit z prawdziwymi nachyleniami, wysokościami i geometrią zderzenia; patrz `docs/VALIDATION.md`.

## Korzystanie bez Unity

**Demo** – co wydarzyło się w scenariuszu i jak rozprzestrzenia się chmura odłamków:
```bash
uv run python -m kessler_sim.demo run iridium_cosmos_2009 --times 0,10m,1h,1d,30d --plot
uv run python -m kessler_sim.demo approaches --duration 10m --top 10     # bliskie przeloty w katalogu
```
`--plot` zapisuje wykres w `results/figures/cloud_<scenariusz>.png`.

**Klient testowy** – sprawdza działający serwer bez Unity (w drugim terminalu):
```bash
uv run python -m kessler_sim.client info     # co serwer wysyła po połączeniu
uv run python -m kessler_sim.client watch    # 5 s działania i liczba aktualizacji
uv run python -m kessler_sim.client demo     # najbliższy przelot -> zderzenie -> skok przed i po
```

**Walidacja** – porównuje symulowane odłamki obu zdarzeń historycznych z prawdziwymi odłamkami, które są dziś nadal na orbicie:
```bash
uv run python -m kessler_sim.validation                              # oba zdarzenia
uv run python -m kessler_sim.validation --scenario fengyun_1c_2007   # jedno z nich
```
Zapisuje wykresy w `results/figures/validation_*.png` i tabele w `results/validation_*.json` / `.md`. Wyniki i ich interpretacja są w `docs/VALIDATION.md`.

## Aktualizacja danych

Katalog to zdjęcie pobrane z CelesTrak (obecnie 2026‑11‑12). Aby użyć nowszych danych, wykonaj kolejno:

| Krok | Polecenie | Wynik |
|---|---|---|
| 1. Pobranie | `uv run python src/data_pipeline/TLE_data_download.py` | Aktywne satelity i cztery grupy odłamków w `data/raw/tle/` (wymaga internetu) |
| 2. Sprawdzenie | `uv run python src/data_pipeline/validate_tle.py` | Raport w `data/metadata/tle_validation_report.json` / `_details.csv` |
| 3. Przetworzenie | `uv run python src/data_pipeline/process_TLE_to_json.py --epoch 2026-11-12T22:50:04+00:00` | Wspólny zestaw danych `data/processed/satellites_state.json` z podsumowaniem; chwila zestawu danych (epoka) to podany czas albo bieżąca chwila, jeśli pominie się `--epoch` |
| 4. Zdjęcie dla Unity | `uv run python src/data_pipeline/export_to_unity.py` | Kopiuje zestaw danych do projektu Unity na potrzeby trybu offline (jeśli poprzedni plik się różni, trafia do kopii w `KesslerSimulation/DatasetBackups/`) |

Następnie uruchom serwer ponownie. Ponowne wykonanie kroku 3 z tym samym `--epoch` odtwarza dokładnie ten sam zestaw danych.

## Testy

```bash
uv run pytest        # ok. 145 testów, ok. 12 s
```
Jeden test jest oznaczony jako oczekiwana porażka: dokumentuje znane ograniczenie modelu rozpadu (rzadkie słabe zderzenia).

## Dokładność i ograniczenia

- Pozycje zgadzają się z modelem SGP4 z dokładnością ok. 10 km na niskiej orbicie przez kilka dni (znacznie gorzej dla kilku bardzo wydłużonych orbit dalekiego kosmosu, które są poza zakresem projektu).
- Zderzenia zachodzą tylko na żądanie (Twoje, scenariusza albo z listy bliskich przelotów); nic nie zderza się samo. Planowana jest długoterminowa kaskada ze statystycznie losowanymi zderzeniami przy prawdziwych rozmiarach obiektów.
- Masy i rozmiary obiektów katalogu są szacunkami zależnymi od typu obiektu.
- W porównaniu z prawdziwymi odłamkami z 2009 i 2007 r.: płaszczyzny orbit się zgadzają, ale odrzuty fragmentów są za silne (chmury rozchodzą się szerzej niż w rzeczywistości), fragmenty nie podlegają oporowi atmosfery, a test na Fengyun‑1C daje za mało fragmentów. Szczegóły w `docs/VALIDATION.md`.

## Starsze narzędzia

Te skrypty z wcześniejszych etapów projektu nadal działają; dwa pierwsze zastępuje pakiet `kessler_sim` opisany wyżej:

| Skrypt | Status |
|---|---|
| `src/collision_detection/close_approach.py` (10-minutowe wyszukiwanie bliskich przelotów, zapisuje `data/processed/close_approaches.json`) | zastąpiony przez `kessler_sim.demo approaches` i panel bliskich przelotów w Unity |
| `src/collision_physics/run_breakups.py` (statystyki rozpadu dla bliskich przelotów) | zastąpiony przez zderzenia w `kessler_sim` (demo, serwer, Unity) |
| `src/orbital_mechanics/orbital_elements.py` (elementy orbitalne i wykres wysokość/nachylenie `results/figures/orbital_population_map.png`) | nie zastąpiony – nadal jedyny wykres przeglądowy całej populacji |

## Dalsza dokumentacja

| Dokument | Zawartość |
|---|---|
| `docs/SIMULATION_ENGINE.md` | Jak działa silnik symulacji, jego model i dokładność |
| `docs/VALIDATION.md` | Porównanie z prawdziwymi odłamkami zdarzeń z 2009 i 2007 r. |
| `../KesslerSimulation/README_PL.md` | Wizualizacja w Unity: sterowanie, panele, obserwacja zderzeń |
| `docs/UNITY_AUDIT_EN.md.old`, `docs/UNITY_AUDIT_PL.md.old` | Nieaktualne: przegląd starego prototypu Unity sprzed podłączenia do tego backendu |

Wersja angielska tego pliku to `README.md`.
