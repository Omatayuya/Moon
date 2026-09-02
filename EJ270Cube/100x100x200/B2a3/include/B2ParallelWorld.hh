#ifndef MyParallelWorld_h
#define MyParallelWorld_h 1

#include "G4VUserParallelWorld.hh"
#include "globals.hh"

class B2aDetectorConstruction;

class MyParallelWorld : public G4VUserParallelWorld {
public:
    // コンストラクタで並行世界に名前を付けます
    MyParallelWorld(G4String worldName);
    MyParallelWorld(G4String worldName, const B2aDetectorConstruction* detector);
    virtual ~MyParallelWorld() = default;

    // ジオメトリ（形と配置）を定義するメソッド
    virtual void Construct() override;
    
    // 検出器（スコアラー）を定義するメソッド
    virtual void ConstructSD() override;

private:
    const B2aDetectorConstruction* fDetector = nullptr;
};

#endif
